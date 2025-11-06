import math

a, b = -1.2, 1.2
E0 = 1e-3

def f(x):
    denominator = math.cbrt(((x + 1.2)**2) * ((x - 12)**2))
    if abs(denominator) < 1e-15:
        return math.sin(x + 1.2) / 1e-15
    return math.sin(x + 1.2) / denominator

def simpson(f, a, b, n):
    if n <= 0: return 0.0
    if n % 2 != 0: n += 1

    result = 0.0
    h = (b - a) / n
    for i in range(n):
        x1 = a + i * h
        x2 = a + (i + 1) * h
        mid = (x1 + x2) / 2.0
        try:
            result += (x2 - x1) / 6.0 * (f(x1) + 4*f(mid) + f(x2))
        except ZeroDivisionError:
            continue
    return result

def find_n_by_runge(f, a, b, n0=2):
    p = 4
    n = max(2, int(n0))
    n = n + 1 if n % 2 == 1 else n

    while n <= 2048:
        try:
            I_n = simpson(f, a, b, n)
            I_2n = simpson(f, a, b, 2 * n)
            err_est = abs(I_2n - I_n) / (2**p - 1)

            if err_est <= E0:
                return 2 * n, I_2n, err_est, I_n, I_2n
        except (ZeroDivisionError, ValueError) as e:
            print(f"Ошибка при n={n}: {e}")
            break

        n *= 2
    return n, I_2n, err_est, I_n, I_2n

try:
    n_used, I_approx, err_est, I_n, I_2n = find_n_by_runge(f, a, b, n0=2)
    residual = abs(I_2n - I_n)
    print(f"Использовано отрезков n = {n_used}")
    print(f"I ≈ {I_approx:.9f}")
    print(f"Оценка ошибки по Рунге: {err_est:.9e}")
    print(f"Невязка: {residual:.9e}")
except Exception as e:
    print(f"Ошибка: {e}")