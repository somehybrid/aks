import aks
import pytest


def is_prime(n: int):
    if n < 2:
        return False

    for i in range(2, n):
        if n % i == 0:
            return False
    return True


DEFINITE_PRIMES = [i for i in range(100) if is_prime(i)]
DEFINITE_COMPOSITES = [i for i in range(100) if not is_prime(i)]


@pytest.mark.parametrize("n", DEFINITE_PRIMES)
def test_prime(n: int):
    assert aks.aks(n)


@pytest.mark.parametrize("n", DEFINITE_COMPOSITES)
def test_not_prime(n: int):
    assert not aks.aks(n)
