from __future__ import annotations
from cmath import exp, pi
import copy
import itertools


def fft(a: list[int]) -> list[complex]:
    data: list[complex] = copy.copy(a)
    omega = exp(2 * pi * 1j / len(data))
    double_stride = 1

    b = copy.copy(data)
    while double_stride < len(data):
        stride = double_stride
        double_stride *= 2

        half = len(data) // double_stride
        power = 1

        for i in range(half):
            for start in range(stride):
                b[start + stride * (2 * i)] = (
                    data[start + stride * i] + data[start + stride * (i + half)]
                )
                b[start + stride * (2 * i + 1)] = (
                    data[start + stride * i] - data[start + stride * (i + half)]
                ) * power
            power *= omega

        data, b = b, data
        omega *= omega

    return data


def pad(data: list[int], n: int) -> list[int]:
    diff_length = n - len(data)
    np2 = 1 << n.bit_length()
    return data + [0] * diff_length + [0] * (np2 - n)


def ifft(data: list[complex]) -> list[complex]:
    data = data.copy()
    double_stride = 1
    b = data.copy()

    omega = exp(2 * pi * 1j / len(data))

    while double_stride < len(data):
        stride = double_stride
        double_stride *= 2

        half = len(data) // double_stride
        power = 1

        for i in range(half):
            for start in range(stride):
                b[start + stride * (2 * i)] = (
                    data[start + stride * i] + data[start + stride * (i + half)]
                ) / 2
                b[start + stride * (2 * i + 1)] = (
                    data[start + stride * i] - data[start + stride * (i + half)]
                ) / (2 * power)
            power *= omega

        data, b = b, data
        omega *= omega

    return data


def to_int(data: list[complex]) -> list[int]:
    return [int(round(x.real)) for x in data]


class Polynomial:
    """
    Implements polynomial logic.
    Don't use this for anything other than implementing AKS.
    """

    def __init__(self, coefficients: list[int]):
        self.coefficients = coefficients
        self.degree = len(coefficients) - 1

    def naive_multiply(self, other: Polynomial, r: int, modulus: int) -> Polynomial:
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

                coefficients[deg] %= modulus

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
            return self.naive_multiply(other, r, modulus)

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

    def fft_mul(self, other: Polynomial, r: int, modulus: int) -> Polynomial:
        """
        Calculates polynomial multiplication using the FFT algorithm
        """
        n = self.degree + other.degree + 1

        fft_self = fft(pad(self.coefficients, n))
        fft_other = fft(pad(other.coefficients, n))

        out = [i * j for i, j in zip(fft_self, fft_other)]
        out = to_int(ifft(out))

        out = Polynomial(out)
        out.strip_leading_zeros()
        out = out.reduce(r)
        out %= modulus

        return out

    def multiply(self, other: Polynomial, r: int, modulus: int) -> Polynomial:
        if self.degree + other.degree <= 8:
            return self.naive_multiply(other, r, modulus)
        elif self.degree + other.degree <= 100:
            return self.karatsuba(other, r, modulus)

        return self.fft_mul(other, r, modulus)

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
                output = output.multiply(base, r, modulus)
            base = base.multiply(base, r, modulus)
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
