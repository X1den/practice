import math

def w_poly(x, nodes):
    result = 1.0
    for xi in nodes:
        result *= (x - xi)
    return result

def w_prime_at_node(nodes, idx):
    xi = nodes[idx]
    prod = 1.0
    for j, xj in enumerate(nodes):
        if j == idx:
            continue
        prod *= (xi - xj)
    return prod

def integrate_trapezoid(func, a, b, n=1000):
    h = (b - a) / n
    s = 0.5 * (func(a) + func(b))
    for i in range(1, n):
        s += func(a + i * h)
    return s * h

def compute_Cs(a, b, p_func, nodes, n_integral=1000):
    Cs = []
    N = len(nodes)
    for i in range(N):
        xi = nodes[i]
        wprime = w_prime_at_node(nodes, i)
        if wprime == 0.0:
            raise ValueError("Повторяющиеся узлы или ошибка: w'(x_i)=0")

        def integrand(x):
            prod = 1.0
            for j, xj in enumerate(nodes):
                if j != i:
                    prod *= (x - xj)
            return p_func(x) * prod / wprime

        Ci = integrate_trapezoid(integrand, a, b, n=n_integral)
        Cs.append(Ci)
    return Cs

def interpolatory_quadrature(a, b, p_func, f_func, nodes, Cs=None):
    if Cs is None:
        Cs = compute_Cs(a, b, p_func, nodes)
    I_approx = 0.0
    for ci, xi in zip(Cs, nodes):
        I_approx += ci * f_func(xi)
    return I_approx

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
f = lambda x: x**3 * math.exp(-0.12 * x**1.5)

nodes = [1.5, 2.25, 3.0, 3.75]

Cs = compute_Cs(a, b, p, nodes)
print("C_i:", Cs)

I_approx = interpolatory_quadrature(a, b, p, f, nodes, Cs)
print("Приближение интеграла I ≈", I_approx)

M_est = 10.0
bound, integral_absw = error_bound_theoretical(a, b, p, nodes, M_est)
print("∫ p(x)|w(x)| dx =", integral_absw)
print("Теоретическая верхняя граница остатка |R| <=", bound)