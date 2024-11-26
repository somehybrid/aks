from __future__ import annotations
from .polynomial import Polynomial
import math


def phi(n: int) -> int:
    """
    Returns the Euler totient function of n.
    """
    totient = n
    for i in range(2, n):
        if n % i == 0:
            totient *= i - 1
            totient //= i
            while n % i == 0:
                n //= i

    return totient


def is_power(n: int, base: int) -> bool:
    for exponent in range(1, n):
        if base**exponent == n:
            return True
        elif base**exponent > n:
            return False

    return False


def check_power(n: int) -> bool:
    """
    Checks if n is a power of a number.
    """
    for b in range(2, n.bit_length() + 1):
        if is_power(n, b):
            return True

    return False


def order(r: int, n: int) -> int:
    """
    Returns the multiplicative order of n modulo r.
    """

    for i in range(1, n):
        if pow(n, i, r) in [0, 1]:
            return i

    assert False, "unreachable"  # unreachable


# bernstein theorem 2.3
def find_r(n: int) -> tuple[int, int]:
    # this section is based on https://github.com/danaj/Math-Prime-Util-GMP/blob/master/aks.c#L298-L354
    log2n = math.log2(n)
    max_k = int(log2n * log2n)

    r = 2
    s = 0

    while True:
        if n <= r:
            break

        r += 1

        gcd = math.gcd(n, r)

        if gcd != 1:
            return (-1, -1)

        v = order(r, n)
        if v >= max_k:
            continue

        t = phi(r)
        q = t
        phiv = q / v

        slim = 20 * (r - 1)

        t = math.comb(q + slim - 1, slim)
        sbin = math.log2(t)

        if sbin < 2 * math.floor(math.sqrt(q)) * log2n:
            continue

        for s in range(2, slim):
            t = math.comb(q + s - 1, s)
            sbin = math.log2(t)
            if sbin < 2 * math.floor(math.sqrt(q)) * log2n:
                continue

            d = 2
            while d < phiv:
                if phiv % d != 0:
                    d += 1
                    continue
                scmp = 2 * d * math.floor(math.sqrt(q / d)) * log2n

                if sbin < scmp:
                    break

                d += 1

            if d >= phiv:
                break

        if s < slim:
            break

    return r, s


def aks(n: int) -> bool:
    """
    Checks if n is a prime number.
    """
    if n < 31:
        return n in [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]

    if check_power(n):
        return False

    r, s = find_r(n)

    if r == s == -1:
        return False

    if math.gcd(n, r) != 1:
        return False

    if n <= r:
        return True

    for a in range(2, r):
        if n % a == 0:
            return False

    for a in range(2, s):
        poly = Polynomial([a, 1])

        if poly.pow(n, r, n) != Polynomial([a % n] + [0] * (n - 1) + [1]).reduce(r):
            return False

    return True
