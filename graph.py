import numpy as np
import matplotlib.pyplot as plt
import math

def f(x):
    return x**3 * math.exp(-0.12 * x**1.5)

def f_fourth_analytic(x):
    return (81*((81*x**5)/(6250000)+(43*x**2)/(2000))*math.exp(-(3*x**(3/2))/(25)) - 
            (81*math.sqrt(x)*(54*x**3+4375)*math.exp(-(3*x**(3/2))/(25)))/(50000))

def numerical_derivative(f, x, n=1, dx=1e-6):
    if n == 0:
        return f(x)
    elif n == 1:
        return (f(x + dx) - f(x - dx)) / (2 * dx)
    elif n == 2:
        return (f(x + dx) - 2*f(x) + f(x - dx)) / (dx**2)
    elif n == 3:
        return (f(x + 2*dx) - 2*f(x + dx) + 2*f(x - dx) - f(x - 2*dx)) / (2 * dx**3)
    elif n == 4:
        return (f(x + 2*dx) - 4*f(x + dx) + 6*f(x) - 4*f(x - dx) + f(x - 2*dx)) / (dx**4)

x_vals = np.linspace(1.5, 3.75, 1000)

f_vals = [f(x) for x in x_vals]
f_prime_vals = [numerical_derivative(f, x, 1) for x in x_vals]
f_double_prime_vals = [numerical_derivative(f, x, 2) for x in x_vals]
f_triple_prime_vals = [numerical_derivative(f, x, 3) for x in x_vals]
f_fourth_numeric = [numerical_derivative(f, x, 4) for x in x_vals]
f_fourth_analytic_vals = [f_fourth_analytic(x) for x in x_vals]

M_1 = max(abs(df) for df in f_prime_vals)
M_2 = max(abs(df) for df in f_double_prime_vals)
M_3 = max(abs(df) for df in f_triple_prime_vals)
M_4_numeric = max(abs(df) for df in f_fourth_numeric)
M_4_analytic = max(abs(df) for df in f_fourth_analytic_vals)

print("МАКСИМУМЫ ПРОИЗВОДНЫХ НА [1.5, 3.75]:")
print(f"|f'(x)|  max = {M_1:.2f}")
print(f"|f''(x)| max = {M_2:.2f}")
print(f"|f'''(x)| max = {M_3:.2f}")
print(f"|f⁽⁴⁾(x)| численно  max = {M_4_numeric:.2f}")
print(f"|f⁽⁴⁾(x)| аналитически max = {M_4_analytic:.2f}")

diff = max(abs(a - n) for a, n in zip(f_fourth_analytic_vals, f_fourth_numeric))
print(f"\nМаксимальное расхождение между аналитической и численной производной: {diff:.2e}")

plt.figure(figsize=(15, 10))

plt.subplot(2, 3, 1)
plt.plot(x_vals, f_vals, 'b-', linewidth=2)
plt.title('f(x) = x³·exp(-0.12·x¹·⁵)')
plt.grid(True)
plt.xlabel('x')
plt.ylabel('f(x)')

plt.subplot(2, 3, 2)
plt.plot(x_vals, f_prime_vals, 'r-', linewidth=2)
plt.title(f"f'(x)\nmax = {M_1:.2f}")
plt.grid(True)
plt.xlabel('x')
plt.ylabel("f'(x)")

plt.subplot(2, 3, 3)
plt.plot(x_vals, f_double_prime_vals, 'g-', linewidth=2)
plt.title(f"f''(x)\nmax = {M_2:.2f}")
plt.grid(True)
plt.xlabel('x')
plt.ylabel("f''(x)")

plt.subplot(2, 3, 4)
plt.plot(x_vals, f_triple_prime_vals, 'm-', linewidth=2)
plt.title(f"f'''(x)\nmax = {M_3:.2f}")
plt.grid(True)
plt.xlabel('x')
plt.ylabel("f'''(x)")

plt.subplot(2, 3, 5)
plt.plot(x_vals, f_fourth_analytic_vals, 'c-', linewidth=2, label='Аналитическая')
plt.plot(x_vals, f_fourth_numeric, 'r--', linewidth=1, label='Численная')
plt.title(f"f⁽⁴⁾(x) аналитически\nmax = {M_4_analytic:.2f}")
plt.grid(True)
plt.xlabel('x')
plt.ylabel("f⁽⁴⁾(x)")
plt.legend()

plt.subplot(2, 3, 6)
plt.plot(x_vals, f_prime_vals, 'r-', label="f'(x)", linewidth=1)
plt.plot(x_vals, f_double_prime_vals, 'g-', label="f''(x)", linewidth=1)
plt.plot(x_vals, f_triple_prime_vals, 'm-', label="f'''(x)", linewidth=1)
plt.plot(x_vals, f_fourth_analytic_vals, 'c-', label="f⁽⁴⁾(x)", linewidth=2)
plt.title('Все производные')
plt.grid(True)
plt.xlabel('x')
plt.ylabel('Производные')
plt.legend()

plt.tight_layout()
plt.show()

def w_poly(x, nodes):
    result = 1.0
    for xi in nodes:
        result *= (x - xi)
    return result

def integrate_trapezoid(func, a, b, n=1000):
    h = (b - a) / n
    s = 0.5 * (func(a) + func(b))
    for i in range(1, n):
        s += func(a + i * h)
    return s * h

def error_bound_theoretical(a, b, p_func, nodes, M_np1, n_integral=1000):
    def integrand(x):
        return abs(p_func(x) * w_poly(x, nodes))
    integral_absw = integrate_trapezoid(integrand, a, b, n=n_integral)
    N = len(nodes) - 1
    bound = M_np1 * integral_absw / math.factorial(N + 1)
    return bound, integral_absw

a = 1.5
b = 3.75
p = lambda x: 1.0
nodes = [1.5, 2.25, 3.0, 3.75]

integral_absw = error_bound_theoretical(a, b, p, nodes, M_4_analytic)[1]
bound_correct = M_4_analytic * integral_absw / math.factorial(4)

print(f"\nОЦЕНКА ПОГРЕШНОСТИ:")
print(f"M₄ (аналитически) = {M_4_analytic:.2f}")
print(f"∫ p(x)|w(x)| dx = {integral_absw:.6f}")
print(f"Правильная граница |R| <= {bound_correct:.6f}")
print(f"Старая граница с M_est=10: {10.0 * integral_absw / math.factorial(4):.6f}")
print(f"Занижение оценки в {M_4_analytic/10:.1f} раз!")