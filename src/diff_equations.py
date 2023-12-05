from typing import Callable

# TODO: POSIBLE  TYPE ALIASING
RealFunction = Callable[[float, ...], float]


# --------------- EULER INITIAL-VALUE PROBLEM ----------#

def euler(f: Callable[[float, float], float],
          a: float, b: float,
          n: int, alpha: float) -> tuple[list, list]:
    t = [0.0] * n
    w = [0.0] * n

    h = (b - a) / n
    t[0] = a
    w[0] = alpha

    for i in range(1, n):
        w[i] = w[i - 1] + h * f(t[i - 1], w[i - 1])
        t[i] = a + i * h

    return t, w


def runge_kutta(f: Callable[[float, float], float],
                a: float, b: float, ):
   pass
