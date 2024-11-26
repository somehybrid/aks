import time
from scipy.interpolate import interp1d
from scipy.signal import savgol_filter
import numpy
import typing
from collections.abc import Iterable
from concurrent.futures import ProcessPoolExecutor
from matplotlib import pyplot as plt
import aks
import math


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

    for i in range(2, int(math.sqrt(n)) + 1):
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


def plot_benches():
    n = 2000
    k1 = []
    v1 = []
    for i in range(1, n):
        _ = is_prime(i)  # warmup
    for i in range(1, n):
        if i < 31:
            _ = is_prime(i)  # warmup
        if is_prime(i):
            k1.append(i)
            v1.append(bench(i, is_prime))
    print("finished naive prime test")

    workers = 50
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

    print("finished aks primality test")
    k1 = numpy.array(k1)
    v1 = numpy.array(v1)
    k2 = numpy.array(k2)
    v2 = numpy.array(v2)

    interpolator = interp1d(k1, v1, kind="cubic")
    interpolator2 = interp1d(k2, v2, kind="cubic")

    x1 = numpy.linspace(k1.min(), k1.max(), num=50)
    x2 = numpy.linspace(k2.min(), k2.max(), num=50)

    y1 = interpolator(x1)
    y2 = interpolator2(x2)

    coefs1 = numpy.polyfit(x1, y1, 4)

    poly1 = numpy.poly1d(coefs1)

    coefs1 = numpy.polyfit(x1, y1, 4)
    poly1 = numpy.poly1d(coefs1)
    coefs2 = numpy.polyfit(x2, y2, 4)
    poly2 = numpy.poly1d(coefs2)

    figure, axis = plt.subplots(1, 2)
    axis[0].set_xscale("log")
    axis[0].plot(x1, y1, label="time")
    axis[0].plot(x1, poly1(x1), label="polynomial fit")
    axis[0].legend()
    axis[0].set_title("naive prime test (testing up to sqrt(n))")

    axis[1].set_xscale("log")
    axis[1].plot(x2, y2, label="time")
    axis[1].plot(x2, poly2(x2), label="polynomial fit")
    axis[1].legend()
    axis[1].set_title("aks primality test")
    plt.show()


if __name__ == "__main__":
    plot_benches()
