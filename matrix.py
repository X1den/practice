import math

a, b = -1.2, 1.2
E0 = 1e-3

def f(x):
    denom = ((x + 1.2)**2) * ((x - 12)**2)
    if abs(denom) < 1e-15:
        return 0.0
    return math.sin(x + 1.2) / math.cbrt(denom)

def simpson(f, a, b, n):
    if n <= 0: return 0.0
    if n % 2 != 0: n += 1
    h = (b - a) / n
    result = 0.0
    for i in range(n):
        x1 = a + i * h
        x2 = a + (i + 1) * h
        mid = (x1 + x2) / 2
        try:
            result += (x2 - x1) / 6 * (f(x1) + 4*f(mid) + f(x2))
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
        except (ZeroDivisionError, ValueError):
            break
        n *= 2
    return n, I_2n, err_est, I_n, I_2n

target = 2 * E0
eps = 1e-6
max_eps = 1e-2

while eps < max_eps:
    try:
        val = abs(simpson(f, a + eps/100, a + eps, 50))
        if val <= target:
            break
    except Exception:
        pass
    eps *= 1.5

a_adj = a + eps

def main():
    n_used, I_approx, err_est, I_n, I_2n = find_n_by_runge(f, a_adj, b, n0=2)
    residual = abs(I_2n - I_n)
    print(f"Найдено eps ≈ {eps:.6e}")
    print(f"Интеграл на [-1.2, -1.2+eps] ≈ {val:.6e} (<= {target})")
    print(f"Использовано отрезков n = {n_used}")
    print(f"I ≈ {I_approx:.9f}")
    print(f"Оценка ошибки по Рунге: {err_est:.9e}")
    print(f"Невязка: {residual:.9e}")
    print(f"Интегрирование выполнено на отрезке [{a_adj:.6f}, {b}] (смещение {a_adj - a:.6f})")

if __name__ == "__main__":
    main()