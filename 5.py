import math

def f(x):
    return x**3 * math.exp(-0.12 * x**1.5)

def P(n, x):
    if n == 0: return 1.0
    if n == 1: return x
    p0, p1 = 1.0, x
    for k in range(1, n):
        p2 = ((2*k + 1)*x*p1 - k*p0) / (k + 1)
        p0, p1 = p1, p2
    return p1

def P_prime(n, x):
    if abs(x*x - 1) < 1e-12: return 0.0
    return n * (x * P(n, x) - P(n-1, x)) / (x*x - 1)

def gauss_legendre(n):
    nodes, weights = [], []
    m = (n + 1) // 2
    for i in range(m):
        x = math.cos(math.pi * (i + 0.75) / (n + 0.5))
        while True:
            pn = P(n, x)
            pn_prime = P_prime(n, x)
            if abs(pn_prime) < 1e-12: break
            dx = pn / pn_prime
            x -= dx
            if abs(dx) < 1e-16: break
        weight = 2.0 / ((1 - x*x) * pn_prime**2)
        nodes.extend([-x, x])
        weights.extend([weight, weight])

    paired = sorted(zip(nodes, weights))
    return [p[0] for p in paired], [p[1] for p in paired]

def w(xi, nodes):
    prod = 1.0
    for xj in nodes:
        prod *= (xi - xj)
    return prod

def int_w2_sq(nodes, N=20000):
    a, b = -1.0, 1.0
    h = (b - a) / N
    integral = 0.0
    for i in range(N):
        x0 = a + i * h
        x1 = x0 + h
        integral += 0.5 * h * (w(x0, nodes)**2 + w(x1, nodes)**2)
    return integral

def quad_gauss(f, a, b, n):
    xi, wi = gauss_legendre(n)
    I = 0.0
    for i in range(n):
        x = 0.5 * (b - a) * xi[i] + 0.5 * (a + b)
        weight = 0.5 * (b - a) * wi[i]
        I += weight * f(x)
    return I

def nth_derivative(f, x, n, h=1e-4):
    if n == 0: return f(x)
    return (nth_derivative(f, x + h, n-1, h) - nth_derivative(f, x - h, n-1, h)) / (2 * h)

def M2n(f, a, b, n_diff, points=100):
    h = (b - a) / points
    max_val = 0.0
    for i in range(points + 1):
        x = a + i * h
        val = nth_derivative(f, x, n_diff, h=1e-5)
        max_val = max(max_val, abs(val))
    return max_val

a, b = 1.5, 3.75
n = 4

xi, wi = gauss_legendre(n)                    # (28) — узлы и веса
I = quad_gauss(f, a, b, n)                     # (28) — квадратура
w2_int = int_w2_sq(xi)                        # ∫_{-1}^1 w²(ξ) dξ
scale = ((b - a) / 2) ** (2*n + 1)             # масштаб для ∫_a^b w²(x) dx
int_w2_ab = w2_int * scale
M8 = M2n(f, a, b, 8, points=50)               # M_{2n}
R_bound = M8 / math.factorial(2*n) * int_w2_ab # (32)

print(f"=== КНИГА: Теорема 1, формулы (28)–(32) ===")
print(f"n = {n}")
print(f"∫_{a}^{b} f(x) dx ≈ {I:.12f}          ← (28)")
print(f"∫_{-1}^1 w²(ξ) dξ = {w2_int:.6e}")
print(f"∫_{a}^{b} w²(x) dx = {int_w2_ab:.6e}   ← масштаб (b-a)/2")
print(f"M_{2*n} = max |f^{({2*n})}(x)| ≈ {M8:.6e}")
print(f"|R| ≤ M_{2*n} / (2n)! * ∫ w² dx = {R_bound:.6e}  ← (32)")