import math
from decimal import Decimal, getcontext
getcontext().prec = 30

a, b = Decimal('1.5'), Decimal('3.75')
E0 = Decimal('1e-8')

def f(x):
    # x - Decimal
    x_float = float(x)                # для возведения в степень и экспоненты
    return Decimal(x_float**3) * Decimal(math.exp(-0.12 * x_float**1.5))

def simpson(f, a, b, n):
    if n <= 0: return Decimal(0)
    if n % 2 != 0: n += 1

    result = Decimal(0)
    h = (b - a) / Decimal(n)
    for i in range(n):
        x1 = a + i * h
        x2 = a + (i + 1) * h
        mid = (x1 + x2) / Decimal(2)
        result += (x2 - x1) / Decimal(6) * (f(x1) + 4*f(mid) + f(x2))
    return result

def find_n_by_runge(f, a, b, n0=2):
    p = 4
    n = max(2, int(n0))
    n = n + 1 if n % 2 == 1 else n

    while n <= 2048:
        I_n = simpson(f, a, b, n)
        I_2n = simpson(f, a, b, 2 * n)
        err_est = abs(I_2n - I_n) / (Decimal(2)**p - 1)

        if err_est <= E0:
            return 2 * n, float(I_2n), float(err_est), float(I_n), float(I_2n)

        n *= 2
    return n, float(I_2n), float(err_est), float(I_n), float(I_2n)

n_used, I_approx, err_est, I_n, I_2n = find_n_by_runge(f, a, b, n0=2)
residual = abs(I_2n - I_n)
print(f"Использовано отрезков n = {n_used}")
print(f"I ≈ {I_approx:.9f}")
print(f"Оценка ошибки по Рунге: {err_est:.9e}")
print(f"Невязка: {residual:.9e}")