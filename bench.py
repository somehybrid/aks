from collections.abc import Iterable
from concurrent.futures import ProcessPoolExecutor
import typing
import time
import math

from scipy.interpolate import interp1d
from matplotlib import pyplot as plt
from matplotlib.axes import Axes
import numpy
import aks


def bench(n: int, func: typing.Callable[[int], bool]) -> int:
    start = time.perf_counter_ns()
    _ = func(n)
    end = time.perf_counter_ns()

    start2 = time.perf_counter_ns()
    _ = func(n)
    end2 = time.perf_counter_ns()

    return ((end - start) + (end2 - start2)) // 2


def is_prime(n: int) -> bool:
    if n <= 1:
        return False

    for i in range(2, n):
        if n % i == 0:
            return False
    return True


def batch(data: Iterable[int], func: typing.Callable[[int], bool]) -> dict[int, int]:
    out = {}
    for i in data:
        if i % 100 == 0:
            print(i)
        if is_prime(i):
            out[i] = bench(i, func)
    return out


def plot_bench(data: list[int], kx: list[int], axis: Axes, test: str):
    data = numpy.array(data)
    k = numpy.array(kx)

    x = numpy.linspace(k.min(), k.max(), num=100)

    interpolator = interp1d(k, data, kind="cubic")
    y = interpolator(x)

    polyfit_coefficients = numpy.polyfit(x, y, 4)  # polynomial fit
    polyfit = numpy.poly1d(polyfit_coefficients)

    exponential_coefficients = numpy.polyfit(x, y, 1)  # O(2^n)
    exponential = numpy.poly1d(exponential_coefficients)

    extrapolated_x = numpy.linspace(k.min(), int(k.max() * 1.25), num=120)

    axis.set_xscale("log")

    axis.plot(x, y, label="time")
    axis.plot(x, polyfit(x), label="polynomial fit")
    axis.plot(extrapolated_x, exponential(extrapolated_x), label=r"$O(2^n)$")

    axis.legend()
    axis.set_xlabel("Number tested")
    axis.set_ylabel("Time (ns)")

    axis.set_title(test)


def plot_benches():
    n = 1000
    k1 = []
    v1 = []
    for i in range(1, n):
        if is_prime(i):
            k1.append(i)
            v1.append(bench(i, is_prime))
    print("finished naive prime test")

    workers = 25
    batches = [
        range(max(i * (n // workers), 31), (i + 1) * (n // workers))
        for i in range(workers)
    ]

    k2 = []
    v2 = []

    with ProcessPoolExecutor() as executor:
        for i in executor.map(batch, batches, [aks.aks] * workers):
            k2.extend(i.keys())
            v2.extend(i.values())

    figure, axis = plt.subplots(1, 2)
    plt.title("Primality test benchmarks")

    plot_bench(v1, k1, axis[0], "naive prime test")
    plot_bench(v2, k2, axis[1], "naive prime test")

    plt.plot()
    plt.show()


if __name__ == "__main__":
    plot_benches()
