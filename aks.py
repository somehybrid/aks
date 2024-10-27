from __future__ import annotations
import math


class Polynomial:
    """
    Implements polynomial logic.
    Don't use this for anything other than implementing AKS.
    """

    def __init__(self, coefficients: tuple[int, ...] = (0,)):
        self.coefficients = coefficients

        self.degree = len(coefficients) - 1

    def multiply(self, other: Polynomial, r: int, modulus: int) -> Polynomial:
        """
        Calculates multiplication of the polynomial over the polynomial ring `S = (Z/nZ)[X]/(X^r - 1)`
        """
        coefficients = [0] * r
        for degree, coefficient in enumerate(self.coefficients):
            if coefficient == 0:
                continue

            for other_degree, other_coefficient in enumerate(other.coefficients):
                """
                Implementation detail:

                Here, you don't need to calculate the polynomial remainder of `self` modulo `X^r - 1`, as you can substitute `X^r` to `1` and reduce the degree modulo `r`.
                """
                deg = (degree + other_degree) % r  # reducing polynomial modulo X^r - 1
                coefficients[deg] += coefficient * other_coefficient

                coefficients[deg] %= modulus  # reducing coefficients modulo n

        return Polynomial(tuple(coefficients))

    def __str__(self):
        return " + ".join(
            filter(
                None,
                (
                    f"{coefficient}x^{degree}" if coefficient != 0 else None
                    for degree, coefficient in enumerate(self.coefficients)
                ),
            )
        )

    def pow(self, exponent: int, r: int) -> Polynomial:
        """
        Calculates the polynomial to the power of `exponent` in the polynomial ring `S = (Z/nZ)[X]/(X^r - 1)`
        """
        modulus = exponent

        if exponent == 0:
            return Polynomial((1,))
        elif exponent == 1:
            return self

        base = self
        output = Polynomial((1,))
        while exponent > 0:  # exponentiation by squaring
            if exponent % 2 == 1:
                output = output.multiply(base, r, modulus)  # o * s mod (x^r - 1, n)

            base = base.multiply(base, r, modulus)  # s * s mod (x^r - 1, n)
            exponent //= 2

        return output

    def reduce(self, degree: int) -> Polynomial:
        """
        Reduces the polynomial modulo X^degree - 1
        """
        coefficients = [0] * degree
        for index, coefficient in enumerate(self.coefficients):
            coefficients[index % degree] += coefficient

        return Polynomial(tuple(coefficients))

    def __eq__(self, other: Polynomial) -> bool:
        return self.coefficients == other.coefficients


def phi(n: int):
    """
    Calculates Euler's totient function (or number of coprimes less than n) of n.
    """
    # is it efficient? no. do i care? also no
    amount = 0
    for k in range(1, n):
        if math.gcd(n, k) == 1:
            amount += 1
    return amount


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
    for b in range(2, int(math.log2(n)) + 1):
        if is_power(n, b):
            return True

    return False


def find_r(n: int) -> int:
    """
    Finds the smallest r such that the multiplicative order of r is more than log2(n)^2.
    """
    maxK = int(math.log2(n) ** 2)

    nextR = True

    r = 1
    while nextR:
        r += 1
        nextR = False
        k = 1
        while k <= maxK and nextR is False:
            k += 1
            nextR = pow(n, k, r) in [1, 0]

    return r


def aks(n: int) -> bool:
    """
    Checks if n is a prime number.
    """
    if check_power(n):
        return False

    r = find_r(n)

    if math.gcd(n, r) != 1:
        return False

    if n <= r:
        return True

    for a in range(2, r):
        if n % a == 0:
            return False

    for a in range(2, int(math.sqrt(phi(r)) * math.log2(n))):
        poly = Polynomial((a, 1))

        if poly.pow(n, r) != Polynomial((a,) + (0,) * (n - 1) + (1,)).reduce(r):
            return False

    return True
