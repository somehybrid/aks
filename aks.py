from __future__ import annotations
import itertools
import math


class Polynomial:
    """
    Implements polynomial logic.
    Don't use this for anything other than implementing AKS.
    """

    def __init__(self, coefficients: list[int]):
        self.coefficients = coefficients
        self.degree = len(coefficients) - 1

    def multiply(self, other: Polynomial, r: int, modulus: int) -> Polynomial:
        """
        Calculates multiplication of the polynomial over the polynomial ring `S = (Z/nZ)[X]/(X^r - 1)`
        """
        coefficients = [0] * (min(self.degree + other.degree, r) + 1)
        for degree, coefficient in enumerate(self.coefficients):
            if coefficient == 0:
                continue

            for other_degree, other_coefficient in enumerate(other.coefficients):
                """
                Implementation detail:

                Here, you don't need to calculate the polynomial remainder of `self` modulo `X^r - 1`, as you can substitute `X^r` to `1` and reduce the degree modulo `r`.
                """
                deg = (degree + other_degree) % r
                coefficients[deg] += coefficient * other_coefficient

                coefficients[degree] %= modulus

        return Polynomial(coefficients)

    def __add__(self, other: Polynomial):
        return Polynomial(
            list(
                map(
                    lambda a: a[0] + a[1],
                    itertools.zip_longest(
                        self.coefficients, other.coefficients, fillvalue=0
                    ),
                )
            )
        )

    def __sub__(self, other: Polynomial):
        return Polynomial(
            list(
                map(
                    lambda a: a[0] - a[1],
                    itertools.zip_longest(
                        self.coefficients, other.coefficients, fillvalue=0
                    ),
                )
            )
        )

    def __mod__(self, other: int) -> Polynomial:
        return Polynomial(list(map(lambda a: a % other, self.coefficients)))

    def increase_degree(self, m: int) -> Polynomial:
        """
        Increases the degree of the polynomial by `m`
        """
        return Polynomial([0] * m + self.coefficients)

    def karatsuba(self, other: Polynomial, r: int, modulus: int) -> Polynomial:
        """
        Calculates the polynomial multiplication using the Karatsuba algorithm
        """
        m = (max(self.degree, other.degree) + 1) // 2

        if m <= 4:
            return self.multiply(other, r, modulus)

        low1 = Polynomial(self.coefficients[:m])
        high1 = Polynomial(self.coefficients[m:])

        low2 = Polynomial(other.coefficients[:m])
        high2 = Polynomial(other.coefficients[m:])

        # stripping the zeros improves performance
        low1.strip_leading_zeros()
        high1.strip_leading_zeros()
        low2.strip_leading_zeros()
        high2.strip_leading_zeros()

        z0 = low1.karatsuba(low2, r, modulus)
        z2 = high1.karatsuba(high2, r, modulus)
        z3 = (low1 + high1).karatsuba(low2 + high2, r, modulus)

        z1 = z3 - z2 - z0

        result = z2.increase_degree(m * 2) + z1.increase_degree(m) + z0
        result = result.reduce(r) % modulus
        result.strip_leading_zeros()

        return result

    def pow(self, exponent: int, r: int, modulus: int) -> Polynomial:
        """
        Calculates the polynomial to the power of `exponent` in the polynomial ring `S = (Z/nZ)[X]/(X^r - 1)`
        """
        if exponent == 0:
            return Polynomial([1])
        elif exponent == 1:
            return self

        output = Polynomial([1])
        base = self
        while exponent > 0:
            if exponent % 2 == 1:
                output = output.karatsuba(base, r, modulus)
            base = base.karatsuba(base, r, modulus)
            exponent //= 2

        return output

    def reduce(self, degree: int) -> Polynomial:
        """
        Reduces the polynomial modulo X^degree - 1
        """
        if self.degree < degree:
            return self

        coefficients = [0] * degree
        for index, coefficient in enumerate(self.coefficients):
            coefficients[index % degree] += coefficient

        return Polynomial(coefficients)

    def strip_leading_zeros(self):
        if self.coefficients == []:
            return

        while self.coefficients[-1] == 0:
            _ = self.coefficients.pop()
            if self.coefficients == []:
                return

        self.degree = len(self.coefficients) - 1

    def __eq__(self, other: Polynomial) -> bool:  # pyright: ignore[reportImplicitOverride, reportIncompatibleMethodOverride]
        self.strip_leading_zeros()
        other.strip_leading_zeros()
        return self.coefficients == other.coefficients

    def __str__(self) -> str:  # pyright: ignore[reportImplicitOverride]
        def format_coefficient(coefficient: int, exponent: int) -> str:
            co = ""

            if exponent == 0 and coefficient != 0:
                co += f"{coefficient}"
            elif exponent != 0 and coefficient != 0:
                co += f"{coefficient}x^{exponent}"
            elif exponent != 0 and coefficient == 1:
                co += f"x^{exponent}"

            return co

        return " + ".join(
            reversed(
                list(
                    format_coefficient(coefficient, exponent)
                    for exponent, coefficient in enumerate(self.coefficients)
                    if coefficient != 0
                )
            )
        )


def phi(n: int):
    """
    Calculates Euler's totient function (or number of coprimes less than n) of n.
    """
    totient = n
    for i in range(2, n):
        if n % i == 0:
            totient *= i - 1
            totient /= i
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
    if n < 31:
        return n in [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]

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
        poly = Polynomial([a, 1])

        if poly.pow(n, r, n) != Polynomial([a % r] + [0] * (n - 1) + [1]).reduce(r):
            return False

    return True
