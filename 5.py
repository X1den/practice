import math
import numpy as np

def f(x):
    return x**3 * math.exp(-0.12 * x**1.5)

#(28)
GAUSS4_X = [-0.8611363115940526,
            -0.3399810435848563,
             0.3399810435848563,
             0.8611363115940526]

GAUSS4_C = [0.3478548451374538,
            0.6521451548625461,
            0.6521451548625461,
            0.3478548451374538]

def gauss4_integrate(g, a, b):
    s = 0.0
    for xi, ci in zip(GAUSS4_X, GAUSS4_C):
        t = a + (b-a)*(xi + 1)/2  # масштабирование [-1,1] -> [a,b]
        s += ci * g(t)
    return s * (b-a)/2            # масштабируем вес

#Т1
def w_poly_at(x, nodes):
    prod = 1.0
    for xi in nodes:
        prod *= (x - xi)
    return prod

#Оценка погрешности (32), интеграл
def integrate_trap(func, a, b, N):
    if N <= 0:
        return 0.0
    h = (b - a) / N
    s = 0.5 * (func(a) + func(b))
    x = a + h
    for k in range(1, N):
        s += func(x)
        x += h
    return s * h

def approx_8th_derivative(g, x, h):
    # Центральная разностная формула (9 точек), приближает 8-ю производную.
    v = (g(x-4*h) - 8.0*g(x-3*h) + 28.0*g(x-2*h) - 56.0*g(x-h)
         + 70.0*g(x) - 56.0*g(x+h) + 28.0*g(x+2*h) - 8.0*g(x+3*h) + g(x+4*h))
    return v / (h**8)

def estimate_M8(g, a, b, grid_points=201):
    xs = [a + (b-a)*i/(grid_points-1) for i in range(grid_points)]
    max_abs = 0.0
    for x in xs:
        # выбираем шаг h, чтобы x ± 4h не выходили за [a,b]
        h = (b - a) / (grid_points - 1)
        max_possible_h = min((x - a) / 4.0 if x - a > 0 else 0.0,
                             (b - x) / 4.0 if b - x > 0 else 0.0)
        if max_possible_h <= 0:
            continue
        if h > max_possible_h:
            h = max_possible_h
        if h <= 0:
            continue
        try:
            val = approx_8th_derivative(g, x, h)
        except Exception:
            continue
        if abs(val) > max_abs:
            max_abs = abs(val)
    return max_abs

def gauss4_with_error(g, a, b):
    I = gauss4_integrate(g, a, b)

    mid = 0.5*(a + b)
    half = 0.5*(b - a)
    nodes = [mid + half * xi for xi in GAUSS4_X]

    # (30)
    # (32)
    def integrand_w2(x):
        w = w_poly_at(x, nodes)
        return w * w
    J = integrate_trap(integrand_w2, a, b, N=2000)

    M8 = estimate_M8(g, a, b, grid_points=201)

    factorial_8 = math.factorial(8)
    bound = (M8 / factorial_8) * J

    return I, nodes, J, M8, bound

if __name__ == "__main__":
    a, b = 1.5, 3.75
    I, nodes, J, M8, bound = gauss4_with_error(f, a, b)
    print("Квадратурная формула Гаусса (n=4):")
    print(f" I ≈ {I:.12f}    ← формула (28)")
    print(" Узлы x_i (формула (29)):")
    for xi in nodes:
        print(f"  {xi:.12f}")
    print(f" ∫_a^b ω^2(x) dx (часть формулы (32)) ≈ {J:.6e}")
    print(f" M8 ≈ max |f^(8)(x)| (из формулы (32)) ≈ {M8:.6e}")
    print(f" Теоретическая верхняя граница ошибки |R| ≤ {bound:.6e}")

def gauss_legendre(n):
    J = [[0.0]*n for _ in range(n)]
    for i in range(n-1):
        a = i + 1
        J[i][i+1] = a / math.sqrt(4*a*a - 1)
        J[i+1][i] = J[i][i+1]

    from scipy.linalg import eigh 

    eigvals, eigvecs = eigh(J)
    nodes = eigvals.tolist()
    weights = [(2 * (eigvecs[0,i]**2)) for i in range(n)]

    return nodes, weights

n = 4
nodes, weights = gauss_legendre(n)

print("Узлы x_i:")
print(nodes)
print("Весa C_i:")
print(weights)