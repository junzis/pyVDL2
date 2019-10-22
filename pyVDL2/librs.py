# Copyright (c) 2012-2015 Tomer Filiba
# Copyright (c) 2015 rotorgit
# Copyright (c) 2015 Stephen Larroque

"""
Reed Solomon
============

A pure-python `universal errors-and-erasures Reed-Solomon Codec:
http://en.wikipedia.org/wiki/Reed%E2%80%93Solomon_error_correction

::

    >>> rs = RSCodec(10)
    >>> rs.encode([1,2,3,4])
    b'\x01\x02\x03\x04,\x9d\x1c+=\xf8h\xfa\x98M'
    >>> rs.encode(b'hello world')
    b'hello world\xed%T\xc4\xfd\xfd\x89\xf3\xa8\xaa'
    >>> rs.decode(b'hello world\xed%T\xc4\xfd\xfd\x89\xf3\xa8\xaa')
    b'hello world'
    >>> rs.decode(b'heXlo worXd\xed%T\xc4\xfdX\x89\xf3\xa8\xaa')     # 3 errors
    b'hello world'
    >>> rs.decode(b'hXXXo worXd\xed%T\xc4\xfdX\x89\xf3\xa8\xaa')     # 5 errors
    b'hello world'
    >>> rs.decode(b'hXXXo worXd\xed%T\xc4\xfdXX\xf3\xa8\xaa')        # 6 errors - fail
    Traceback (most recent call last):
      ...
    ReedSolomonError: Could not locate error

    >>> rs = RSCodec(12)
    >>> rs.encode(b'hello world')
    b'hello world?Ay\xb2\xbc\xdc\x01q\xb9\xe3\xe2='
    >>> rs.decode(b'hello worXXXXy\xb2XX\x01q\xb9\xe3\xe2=')         # 6 errors - ok
    b'hello world'

    You can skip the API and directly use the library as-is. Here's how:

    First you need to init the precomputed tables:
    >> init_tables(0x11d)
    Pro tip: if you get the error: ValueError: byte must be in range(0, 256),
    please check that your prime polynomial is correct for your field.

    Then to encode:
    >> mesecc = rs_encode_msg(mes, n-k)

    To decode:
    >> mes, ecc = rs_correct_msg(mes + ecc, n-k, erase_pos=erase_pos)

    If the decoding fails, it will normally automatically check and raise a
    ReedSolomonError exception that you can handle.

    If you want to manually check if the repaired message is correct:
    >> rsman.check(rmes + recc, k=k)

    Read the sourcecode's comments for more infos about how it works, and for the various parameters you can setup if
    you need to interface with other RS codecs.

"""

import itertools


class ReedSolomonError(Exception):
    pass


# For efficiency, gf_exp[] has size 2 * GF_SIZE
# so that a simple multiplication of two numbers can be resolved without calling % 255.

gf_exp = bytearray([1] * 512)
gf_log = bytearray(256)
field_charac = int(2 ** 8 - 1)


################### GALOIS FIELD ELEMENTS MATHS ###################


def rwh_primes1(n):
    """ Returns  a list of primes < n """
    sieve = [True] * (n / 2)
    for i in range(3, int(n ** 0.5) + 1, 2):
        if sieve[i / 2]:
            sieve[i * i / 2 :: i] = [False] * ((n - i * i - 1) / (2 * i) + 1)
    return [2] + [2 * i + 1 for i in range(1, n / 2) if sieve[i]]


def find_prime_polys(generator=2, c_exp=8, fast_primes=False, single=False):
    """Compute the list of prime polynomials for the given generator
    and galois field characteristic exponent."""

    root_charac = 2  # GF(2)
    field_charac = int(root_charac ** c_exp - 1)
    field_charac_next = int(root_charac ** (c_exp + 1) - 1)

    prim_candidates = []
    if fast_primes:
        prim_candidates = rwh_primes1(field_charac_next)
        prim_candidates = [x for x in prim_candidates if x > field_charac]
    else:
        prim_candidates = range(
            field_charac + 2, field_charac_next, root_charac
        )

    correct_primes = []

    for prim in prim_candidates:
        seen = bytearray(field_charac + 1)
        conflict = False

        x = 1
        for i in range(field_charac):
            # Compute the next value in the field
            x = gf_mult_noLUT(x, generator, prim, field_charac + 1)

            # Rejection criterion
            if x > field_charac or seen[x] == 1:
                conflict = True
                break

            # Else flag this value as seen
            else:
                seen[x] = 1

        if not conflict:
            correct_primes.append(prim)
            if single:
                return prim

    return correct_primes


def init_tables(prim=0x11D, generator=2, c_exp=8):
    """Precompute the logarithm and anti-log tables for faster computation
    later, using the provided primitive polynomial.These tables are used
    for multiplication/division since addition/substraction are simple XOR
    operations inside GF of characteristic 2.

    More info:
    https://en.wikipedia.org/wiki/Finite_field_arithmetic#Implementation_tricks
    """
    global gf_exp, gf_log, field_charac

    field_charac = int(2 ** c_exp - 1)

    gf_exp = bytearray(field_charac * 2)
    gf_log = bytearray(field_charac + 1)

    x = 1

    for i in range(field_charac):
        gf_exp[i] = x
        gf_log[x] = i

        if generator == 2:  # faster
            x <<= 1  # multiply by 2
            if x & 0x100:  # x >= 256
                # substract the primary polynomial to the current value
                x ^= prim

        else:
            x = gf_mult_noLUT(x, generator, prim, field_charac + 1)

    for i in range(field_charac, field_charac * 2):
        gf_exp[i] = gf_exp[i - field_charac]

    return gf_log, gf_exp


def gf_add(x, y):
    return x ^ y


def gf_sub(x, y):
    # the same as addition ( mod 2)
    return x ^ y


def gf_neg(x):
    return x


def gf_inverse(x):
    return gf_exp[field_charac - gf_log[x]]


def gf_mul(x, y):
    if x == 0 or y == 0:
        return 0
    return gf_exp[(gf_log[x] + gf_log[y]) % field_charac]


def gf_div(x, y):
    if y == 0:
        raise ZeroDivisionError()
    if x == 0:
        return 0
    return gf_exp[(gf_log[x] + field_charac - gf_log[y]) % field_charac]


def gf_pow(x, power):
    return gf_exp[(gf_log[x] * power) % field_charac]


def gf_mult_noLUT(x, y, prim=0, field_charac_full=256, carryless=True):
    """Galois Field integer multiplication using Russian Peasant Multiplication algorithm (faster than the standard multiplication + modular reduction).
    If prim is 0 and carryless=False, then the function produces the result for a standard integers multiplication (no carry-less arithmetics nor modular reduction)."""
    r = 0
    while y > 0:
        if y & 1:
            r = r ^ x if carryless else r + x
        y = y >> 1  # equivalent to y // 2
        x = x << 1  # equivalent to x*2
        if prim > 0 and x & field_charac_full:
            x = x ^ prim

    return r


################### GALOIS FIELD POLYNOMIALS MATHS ###################


def gf_poly_scale(p, x):
    return bytearray([gf_mul(p[i], x) for i in range(len(p))])


def gf_poly_add(p, q):
    r = bytearray(max(len(p), len(q)))
    r[len(r) - len(p) : len(r)] = p

    for i in range(len(q)):
        r[i + len(r) - len(q)] ^= q[i]
    return r


def gf_poly_mul(p, q):
    """Multiply two polynomials, inside Galois Field (but the procedure is generic). Optimized function by precomputation of log."""
    r = bytearray(len(p) + len(q) - 1)

    lp = [gf_log[p[i]] for i in range(len(p))]

    for j in range(len(q)):
        qj = q[j]
        if qj != 0:
            lq = gf_log[qj]
            for i in range(len(p)):
                if p[i] != 0:
                    r[i + j] ^= gf_exp[lp[i] + lq]
    return r


def gf_poly_neg(poly):
    """Returns the polynomial with all coefficients negated.

    In GF(2^p), negation does not change the coefficient,  return the polynomial as-is.
    """
    return poly


def gf_poly_div(dividend, divisor):
    """Fast polynomial division.

    Using Extended Synthetic Division and optimized for GF(2^p) computations
    (doesn't work with standard polynomials outside of this galois field).

    eg: 1 + 2x + 5x^2 = [5, 2, 1], NOT [1, 2, 5]
    """

    # Copy the dividend list and pad with 0 where the ecc bytes will be computed
    msg_out = bytearray(dividend)

    # normalizer = divisor[0] # precomputing for performance
    for i in range(len(dividend) - (len(divisor) - 1)):
        coef = msg_out[i]
        if coef != 0:
            for j in range(1, len(divisor)):
                if divisor[j] != 0:
                    msg_out[i + j] ^= gf_mul(divisor[j], coef)

    separator = -(len(divisor) - 1)

    return msg_out[:separator], msg_out[separator:]


def gf_poly_square(poly):
    """Linear time implementation of polynomial squaring.

    See paper for detail:
    A fast software implementation for arithmetic operations in GF (2n)
    """
    length = len(poly)
    out = bytearray(2 * length - 1)
    for i in range(length - 1):
        p = poly[i]
        k = 2 * i
        if p != 0:
            out[k] = gf_exp[2 * gf_log[p]]
    out[2 * length - 2] = gf_exp[2 * gf_log[poly[length - 1]]]
    if out[0] == 0:
        out[0] = 2 * poly[1] - 1
    return out


def gf_poly_eval(poly, x):
    """Evaluates a polynomial in GF(2^p) given the value for x.
    This is based on Horner's scheme for maximum efficiency."""
    y = poly[0]
    for i in range(1, len(poly)):
        y = gf_mul(y, x) ^ poly[i]
    return y


################### REED-SOLOMON ENCODING ###################


def rs_generator_poly(nsym, fcr=0, generator=2):
    """Generate an irreducible generator polynomial
    (necessary to encode a message into Reed-Solomon)"""
    g = bytearray([1])
    for i in range(nsym):
        g = gf_poly_mul(g, [1, gf_pow(generator, i + fcr)])
    return g


def rs_generator_poly_all(max_nsym, fcr=0, generator=2):
    """Generate all irreducible generator polynomials up to max_nsym
    (usually you can use n, the length of the message+ecc).
    Very useful to reduce processing time if you want to encode using
    variable schemes and nsym rates.
    """
    g_all = {}
    g_all[0] = g_all[1] = [1]
    for nsym in range(max_nsym):
        g_all[nsym] = rs_generator_poly(nsym, fcr, generator)
    return g_all


def rs_encode_msg(msg_in, nsym, fcr=0, generator=2, gen=None):
    """Reed-Solomon main encoding function, using polynomial division."""
    global field_charac

    if (len(msg_in) + nsym) > field_charac:
        raise ValueError(
            "Message is too long (%i when max is %i)"
            % (len(msg_in) + nsym, field_charac)
        )
    if gen is None:
        gen = rs_generator_poly(nsym, fcr, generator)

    msg_in = bytearray(msg_in)
    msg_out = bytearray(msg_in) + bytearray(len(gen) - 1)

    # Precompute the logarithm of every items in the generator
    lgen = bytearray([gf_log[gen[j]] for j in range(len(gen))])

    # Extended synthetic division main loop
    # Fastest implementation with PyPy (but the Cython version in creedsolo.pyx is about 2x faster)

    for i in range(len(msg_in)):
        coef = msg_out[i]
        if coef != 0:
            lcoef = gf_log[coef]

            for j in range(1, len(gen)):
                msg_out[i + j] ^= gf_exp[lcoef + lgen[j]]

    msg_out[: len(msg_in)] = msg_in
    return msg_out


################### REED-SOLOMON DECODING ###################


def rs_calc_syndromes(msg, nsym, fcr=0, generator=2):
    """Computes the syndromes polynomial.

    Mathematically, it's essentially equivalent to a Fourrier Transform
    (Chien search being the inverse).
    """
    # add a 0 coefficient for the lowest degree
    syndromes = [0] + [
        gf_poly_eval(msg, gf_pow(generator, i + fcr)) for i in range(nsym)
    ]

    return syndromes


def rs_correct_errata(msg_in, synd, err_pos, fcr=0, generator=2):
    """Forney algorithm.

    Computes the values (error magnitude) to correct the input message.
    err_pos is a list of the positions of the errors/erasures/errata
    """

    global field_charac

    msg = bytearray(msg_in)

    coef_pos = [len(msg) - 1 - p for p in err_pos]

    err_loc = rs_find_errata_locator(coef_pos, generator)

    err_eval = rs_find_error_evaluator(synd[::-1], err_loc, len(err_loc) - 1)[
        ::-1
    ]

    X = []
    for i in range(len(coef_pos)):
        l = field_charac - coef_pos[i]
        X.append(gf_pow(generator, -l))

    # Forney algorithm: compute the magnitudes

    E = bytearray(len(msg))
    Xlength = len(X)
    for i, Xi in enumerate(X):

        Xi_inv = gf_inverse(Xi)

        err_loc_prime_tmp = []
        for j in range(Xlength):
            if j != i:
                err_loc_prime_tmp.append(gf_sub(1, gf_mul(Xi_inv, X[j])))

        err_loc_prime = 1
        for coef in err_loc_prime_tmp:
            err_loc_prime = gf_mul(err_loc_prime, coef)

        y = gf_poly_eval(err_eval[::-1], Xi_inv)
        y = gf_mul(gf_pow(Xi, 1 - fcr), y)

        magnitude = gf_div(y, err_loc_prime)

        E[err_pos[i]] = magnitude

    # Apply the correction of values to get our message corrected
    msg = gf_poly_add(msg, E)

    return msg


def rs_find_error_locator(synd, nsym, erase_loc=None, erase_count=0):
    """Find error/errata locator and evaluator polynomials with Berlekamp-Massey algorithm"""

    # Init the polynomials
    if erase_loc:
        err_loc = bytearray(erase_loc)
        old_loc = bytearray(erase_loc)
    else:
        err_loc = bytearray([1])
        old_loc = bytearray([1])

    synd_shift = 0
    if len(synd) > nsym:
        synd_shift = len(synd) - nsym

    for i in range(nsym - erase_count):
        if erase_loc:
            K = erase_count + i + synd_shift
        else:
            K = i + synd_shift

        delta = synd[K]
        for j in range(1, len(err_loc)):
            delta ^= gf_mul(err_loc[-(j + 1)], synd[K - j])

        # Shift polynomials to compute the next degree
        old_loc = old_loc + bytearray([0])

        # Iteratively estimate the errata locator and evaluator polynomials
        if delta != 0:
            if len(old_loc) > len(err_loc):
                new_loc = gf_poly_scale(old_loc, delta)
                old_loc = gf_poly_scale(err_loc, gf_inverse(delta))
                err_loc = new_loc

            # Update with the discrepancy
            err_loc = gf_poly_add(err_loc, gf_poly_scale(old_loc, delta))

    # Check if the result is correct, not too many errors to correct

    # drop leading 0s, else errs will not be of the correct size
    err_loc = list(itertools.dropwhile(lambda x: x == 0, err_loc))
    errs = len(err_loc) - 1
    if (errs - erase_count) * 2 + erase_count > nsym:
        raise ReedSolomonError("Too many errors to correct")

    return err_loc


def rs_find_errata_locator(e_pos, generator=2):
    """Compute the erasures/errors/errata locator polynomial."""
    e_loc = [1]
    for i in e_pos:
        e_loc = gf_poly_mul(e_loc, gf_poly_add([1], [gf_pow(generator, i), 0]))
    return e_loc


def rs_find_error_evaluator(synd, err_loc, nsym):
    """Compute the error."""
    # Omega(x) = [ Synd(x) * Error_loc(x) ] mod x^(n-k+1)
    _, remainder = gf_poly_div(
        gf_poly_mul(synd, err_loc), ([1] + [0] * (nsym + 1))
    )

    return remainder


def rs_find_errors(err_loc, nmess, generator=2):
    """Find the roots (ie, where evaluation = zero) of error polynomial by bruteforce trial, this is a sort of Chien's search (but less efficient, Chien's search is a way to evaluate the polynomial such that each evaluation only takes constant time)."""
    # nmess = length of whole codeword (message + ecc symbols)
    errs = len(err_loc) - 1
    err_pos = []
    for i in range(nmess):
        if gf_poly_eval(err_loc, gf_pow(generator, i)) == 0:
            err_pos.append(nmess - 1 - i)

    # Sanity check: the number of errors/errata positions found should be exactly the same as the length of the errata locator polynomial

    if len(err_pos) != errs:
        raise ReedSolomonError(
            "Too many (or few) errors found by Chien Search for the errata locator polynomial!"
        )
    return err_pos


def rs_forney_syndromes(synd, pos, nmess, generator=2):
    """Compute Forney syndromes.

    It computes a modified syndromes to compute only errors (erasures are trimmed out)
    """
    erase_pos_reversed = [nmess - 1 - p for p in pos]

    # Optimized method, all operations are inlined
    fsynd = list(synd[1:])
    for i in range(len(pos)):
        x = gf_pow(generator, erase_pos_reversed[i])
        for j in range(len(fsynd) - 1):
            fsynd[j] = gf_mul(fsynd[j], x) ^ fsynd[j + 1]

    return fsynd


def rs_correct_msg(
    msg_in, nsym, fcr=0, generator=2, erase_pos=None, only_erasures=False
):
    """Reed-Solomon main decoding function"""
    global field_charac
    if len(msg_in) > field_charac:
        raise ValueError(
            "Message is too long (%i when max is %i)"
            % (len(msg_in), field_charac)
        )

    msg_out = bytearray(msg_in)  # make a copy

    if erase_pos is None:
        erase_pos = []
    else:
        for e_pos in erase_pos:
            msg_out[e_pos] = 0

    if len(erase_pos) > nsym:
        raise ReedSolomonError("Too many erasures to correct")

    synd = rs_calc_syndromes(msg_out, nsym, fcr, generator)

    if max(synd) == 0:
        return msg_out[:-nsym], msg_out[-nsym:]  # no errors

    if only_erasures:
        err_pos = []
    else:
        # compute the Forney syndromes
        fsynd = rs_forney_syndromes(synd, erase_pos, len(msg_out), generator)

        # compute the error locator polynomial using Berlekamp-Massey
        err_loc = rs_find_error_locator(
            fsynd, nsym, erase_count=len(erase_pos)
        )

        # locate the message errors using Chien search (or bruteforce search)
        err_pos = rs_find_errors(err_loc[::-1], len(msg_out), generator)

        if err_pos is None:
            raise ReedSolomonError("Could not locate error")

    # Find errors values and apply them to correct the message
    msg_out = rs_correct_errata(
        msg_out, synd, (erase_pos + err_pos), fcr, generator
    )

    # check if the final message is fully repaired
    synd = rs_calc_syndromes(msg_out, nsym, fcr, generator)

    if max(synd) > 0:
        raise ReedSolomonError("Could not correct message")

    return (msg_out[:-nsym], msg_out[-nsym:])


def rs_correct_msg_nofsynd(
    msg_in, nsym, fcr=0, generator=2, erase_pos=None, only_erasures=False
):
    """Reed-Solomon main decoding function, without using the modified Forney syndromes"""
    global field_charac
    if len(msg_in) > field_charac:
        raise ValueError(
            "Message is too long (%i when max is %i)"
            % (len(msg_in), field_charac)
        )

    msg_out = bytearray(msg_in)  # copy of message

    if erase_pos is None:
        erase_pos = []
    else:
        for e_pos in erase_pos:
            msg_out[e_pos] = 0

    if len(erase_pos) > nsym:
        raise ReedSolomonError("Too many erasures to correct")

    synd = rs_calc_syndromes(msg_out, nsym, fcr, generator)

    if max(synd) == 0:
        return msg_out[:-nsym], msg_out[-nsym:]  # no errors

    erase_loc = None

    erase_count = 0
    if erase_pos:
        erase_count = len(erase_pos)
        erase_pos_reversed = [len(msg_out) - 1 - eras for eras in erase_pos]
        erase_loc = rs_find_errata_locator(
            erase_pos_reversed, generator=generator
        )

    if only_erasures:
        err_loc = erase_loc[::-1]
    else:
        err_loc = rs_find_error_locator(
            synd, nsym, erase_loc=erase_loc, erase_count=erase_count
        )
        err_loc = err_loc[::-1]

    err_pos = rs_find_errors(err_loc, len(msg_out), generator)

    if err_pos is None:
        raise ReedSolomonError("Could not locate error")

    msg_out = rs_correct_errata(
        msg_out, synd, err_pos, fcr=fcr, generator=generator
    )

    synd = rs_calc_syndromes(msg_out, nsym, fcr, generator)
    if max(synd) > 0:
        raise ReedSolomonError("Could not correct message")

    return msg_out[:-nsym], msg_out[-nsym:]


def rs_check(msg, nsym, fcr=0, generator=2):
    # Returns true if the message + ecc has no error of false otherwise
    return max(rs_calc_syndromes(msg, nsym, fcr, generator)) == 0


class RSCodec(object):
    """
    A Reed Solomon encoder/decoder. After initializing the object,
    use ``encode`` to encode a (byte)string to include the RS correction code,
    and pass such an encoded (byte)string to``decode`` to extract the
    original message (if the number of errors allows for correct decoding).

    The ``nsym`` argument is the length of the correction code, and it
    determines the number of error bytes (if I understand this correctly,
    half of ``nsym`` is correctable)
    """

    # Modifications by rotorgit 2/3/2015:
    # Added support for US FAA ADSB UAT RS FEC, by allowing user to specify
    # different primitive polynomial and non-zero first consecutive root (fcr).
    # For UAT/ADSB use, set fcr=120 and prim=0x187 when instantiating
    # the class; leaving them out will default for previous values (0 and
    # 0x11d)

    def __init__(self, **kwargs):
        # nsym=10, nsize=255, fcr=0, prim=0x11D, generator=2, c_exp=8
        # ):
        """Initialize the Reed-Solomon codec."""
        self.nsym = kwargs.get("nsym", 10)
        self.nsize = kwargs.get("nsize", 255)
        self.fcr = kwargs.get("fcr", 0)  # first consecutive root
        self.prim = kwargs.get("prim", 0x11D)  # prime irreducible polynomial
        self.generator = kwargs.get("generator", 2)
        self.c_exp = kwargs.get("c_exp", 8)

        # Initialize the look-up tables for quick multiplication/division
        init_tables(self.prim, self.generator, self.c_exp)

    def encode(self, data):
        """Encode a message (ie, add the ecc symbols) using Reed-Solomon."""
        if isinstance(data, str):
            data = bytearray(data, "latin-1")
        chunk_size = self.nsize - self.nsym
        enc = bytearray()
        for i in range(0, len(data), chunk_size):
            chunk = data[i : i + chunk_size]
            enc.extend(
                rs_encode_msg(
                    chunk, self.nsym, fcr=self.fcr, generator=self.generator
                )
            )
        return enc

    def decode(self, data, erase_pos=None, only_erasures=False):
        """Repair a message, whatever its size is, by using chunking"""

        if isinstance(data, str):
            data = bytearray(data, "utf-8")

        dec = bytearray()

        for i in range(0, len(data), self.nsize):
            # Split the long message in a chunk
            chunk = data[i : i + self.nsize]

            # Extract the erasures for this chunk
            e_pos = []
            if erase_pos:
                e_pos = [x for x in erase_pos if x <= self.nsize]
                erase_pos = [
                    x - (self.nsize + 1) for x in erase_pos if x > self.nsize
                ]

            # Decode/repair this chunk
            dec.extend(
                rs_correct_msg_nofsynd(
                    chunk,
                    self.nsym,
                    fcr=self.fcr,
                    generator=self.generator,
                    erase_pos=e_pos,
                    only_erasures=only_erasures,
                )[0]
            )

        return dec
