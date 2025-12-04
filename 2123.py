import math

a, b = 1, 2
y0 = 1
n = 10
h = (b - a) / n
epsilon = 1e-10

def f(x, y): return (x**2 + y) / x

def exact_solution(x): return x**2

def trapezoid_integral(f_values, x_values):
    integral = 0
    for i in range(len(x_values) - 1):
        integral += (f_values[i] + f_values[i+1]) * (x_values[i+1] - x_values[i]) / 2
    return integral

# 1. МЕТОД ПОСЛЕДОВАТЕЛЬНЫХ ПРИБЛИЖЕНИЙ (ПИКАРА)
def picard_method():
    x_points = [a + i * h for i in range(n + 1)]
    y_picard = [0] * (n + 1)
    y_picard[0] = y0
    
    def y0_func(x):
        return y0
    
    def y1_func(x):
        if x == a:
            return y0
        x_integral = [a + i * (x - a) / 4 for i in range(5)]
        f_values = [f(xi, y0_func(xi)) for xi in x_integral]
        return y0 + trapezoid_integral(f_values, x_integral)
    
    def y2_func(x):
        if x == a:
            return y0
        x_integral = [a + i * (x - a) / 4 for i in range(5)]
        f_values = [f(xi, y1_func(xi)) for xi in x_integral]
        return y0 + trapezoid_integral(f_values, x_integral)
    
    def y3_func(x):
        if x == a:
            return y0
        x_integral = [a + i * (x - a) / 4 for i in range(5)]
        f_values = [f(xi, y2_func(xi)) for xi in x_integral]
        return y0 + trapezoid_integral(f_values, x_integral)
    
    for i, x in enumerate(x_points):
        y_picard[i] = y3_func(x)
    
    return x_points, y_picard

# 2. МЕТОД ЭЙЛЕРА
def euler_method():
    x_points = [a + i * h for i in range(n + 1)]
    y_euler = [0] * (n + 1)
    y_euler[0] = y0
    
    for i in range(n):
        y_euler[i + 1] = y_euler[i] + h * f(x_points[i], y_euler[i])
    
    return x_points, y_euler

# 3. МЕТОД ТРАПЕЦИЙ ЭЙЛЕРА
def euler_trapezoid_method():
    x_points = [a + i * h for i in range(n + 1)]
    y_trapezoid = [0] * (n + 1)
    y_trapezoid[0] = y0
    
    for i in range(n):
        y_pred = y_trapezoid[i] + h * f(x_points[i], y_trapezoid[i])
        y_trapezoid[i + 1] = y_trapezoid[i] + h * (f(x_points[i], y_trapezoid[i]) + f(x_points[i + 1], y_pred)) / 2
    
    return x_points, y_trapezoid

# 4. МЕТОД РУНГЕ-КУТТЫ 4-ГО ПОРЯДКА
def runge_kutta_4th_order():
    x_points = [a + i * h for i in range(n + 1)]
    y_rk4 = [0] * (n + 1)
    y_rk4[0] = y0
    
    for i in range(n):
        k1 = f(x_points[i], y_rk4[i])
        k2 = f(x_points[i] + h/2, y_rk4[i] + h/2 * k1)
        k3 = f(x_points[i] + h/2, y_rk4[i] + h/2 * k2)
        k4 = f(x_points[i] + h, y_rk4[i] + h * k3)
        y_rk4[i + 1] = y_rk4[i] + h/6 * (k1 + 2*k2 + 2*k3 + k4)
    
    return x_points, y_rk4

# 5. МЕТОД АДАМСА 4-ГО ПОРЯДКА
def adams_method():
    y_adams = [0] * (n + 1)
    x_points, y_adams = runge_kutta_4th_order()
    
    for i in range(3, n):
        f_i = f(x_points[i], y_adams[i])
        f_i_1 = f(x_points[i-1], y_adams[i-1])
        f_i_2 = f(x_points[i-2], y_adams[i-2])
        f_i_3 = f(x_points[i-3], y_adams[i-3])
        y_adams[i+1] = y_adams[i] + h/24 * (55*f_i - 59*f_i_1 + 37*f_i_2 - 9*f_i_3)
    
    return x_points, y_adams

# 6. МЕТОД РУНГЕ-РОМБЕРГА для вычисления y(b) с заданной точностью
def runge_romberg_method():
    current_n = n
    y_prev = 0
    y_current = 0
    error = float('inf')
    results = []
    
    print("\nМетод Рунге-Ромберга:")
    print("n\t\ty(b)\t\tПогрешность")
    
    while error > epsilon:
        h_current = (b - a) / current_n
        x_points = [a + i * h_current for i in range(current_n + 1)]
        y = [0] * (current_n + 1)
        y[0] = y0
        
        for i in range(current_n):
            k1 = f(x_points[i], y[i])
            k2 = f(x_points[i] + h_current/2, y[i] + h_current/2 * k1)
            k3 = f(x_points[i] + h_current/2, y[i] + h_current/2 * k2)
            k4 = f(x_points[i] + h_current, y[i] + h_current * k3)
            y[i + 1] = y[i] + h_current/6 * (k1 + 2*k2 + 2*k3 + k4)
        
        y_current = y[-1]
        results.append((current_n, y_current))
        
        if len(results) >= 2:
            y_h = y_current
            y_2h = results[-2][1]
            error = abs(y_h - y_2h) / 15
        
        print(f"{current_n}\t\t{y_current:.8f}\t{error:.2e}")
        
        if error > epsilon:
            current_n += 1
            y_prev = y_current
    
    if len(results) >= 2:
        y_final = y_current + (y_current - y_prev) / 15
    
    print(f"\nТребуемая точность достигнута при n = {current_n}")
    print(f"y({b}) = {y_current:.8f}")
    print(f"Точное значение: {exact_solution(b):.8f}")
    print(f"Абсолютная погрешность: {abs(y_current - exact_solution(b)):.2e}")
    print()
    
    return y_current, current_n, error

x_picard, y_picard = picard_method()
x_euler, y_euler = euler_method()
x_trapezoid, y_trapezoid = euler_trapezoid_method()
x_rk4, y_rk4 = runge_kutta_4th_order()
x_adams, y_adams = adams_method()

y_rr, n_rr, error_rr = runge_romberg_method()

methods = [
    ("Точное решение", exact_solution(b), 0, 0),
    ("Метод Пикара", y_picard[-1], abs(y_picard[-1] - exact_solution(b)), abs(y_picard[-1] - exact_solution(b))/abs(exact_solution(b))),
    ("Метод Эйлера", y_euler[-1], abs(y_euler[-1] - exact_solution(b)), abs(y_euler[-1] - exact_solution(b))/abs(exact_solution(b))),
    ("Трапеций Эйлера", y_trapezoid[-1], abs(y_trapezoid[-1] - exact_solution(b)), abs(y_trapezoid[-1] - exact_solution(b))/abs(exact_solution(b))),
    ("Рунге-Кутта 4-го порядка", y_rk4[-1], abs(y_rk4[-1] - exact_solution(b)), abs(y_rk4[-1] - exact_solution(b))/abs(exact_solution(b))),
    ("Метод Адамса", y_adams[-1], abs(y_adams[-1] - exact_solution(b)), abs(y_adams[-1] - exact_solution(b))/abs(exact_solution(b))),
    ("Рунге-Ромберг", y_rr, abs(y_rr - exact_solution(b)), abs(y_rr - exact_solution(b))/abs(exact_solution(b))),
]

for name, value, abs_err, rel_err in methods:
    print(f"{name:<25} {value:<15.8f} {abs_err:<15.2e} {rel_err:<15.2e}")