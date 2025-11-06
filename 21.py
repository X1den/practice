import math

def f(x, y):
    return (x**2 + y) / x

def picard_iteration(prev_y_func, x0, y0, x, n):
    if n == 0:
        return y0
    
    def integrand(t):
        return f(t, prev_y_func(t))
    
    return y0 + integrate_trapezoid(integrand, x0, x, 1000)

def integrate_trapezoid(func, a, b, n=1000):
    if a == b:
        return 0.0
    h = (b - a) / n
    s = 0.5 * (func(a) + func(b))
    for i in range(1, n):
        s += func(a + i * h)
    return s * h

def solve_picard(x0, y0, x_points, max_iterations=5):
    results = []
    
    def y0_func(x):
        return y0
    
    current_approx = y0_func
    
    print("Метод последовательных приближений Пикара:")
    print(f"Дифференциальное уравнение: y' = (x² + y)/x")
    print(f"Начальное условие: y({x0}) = {y0}")
    print()
    
    for iteration in range(max_iterations + 1):
        y_values = []
        for x in x_points:
            if iteration == 0:
                y_val = y0
            else:
                y_val = y0 + integrate_trapezoid(lambda t: f(t, current_approx(t)), x0, x, 1000)
            y_values.append(y_val)
        
        print(f"Приближение {iteration}:")
        for i, x in enumerate(x_points):
            print(f"  y({x:.1f}) ≈ {y_values[i]:.6f}")
        print()
        
        results.append(y_values)
        
        if iteration < max_iterations:
            def new_approx(x_val, y_vals=y_values, x_pts=x_points):
                for i in range(len(x_pts) - 1):
                    if x_pts[i] <= x_val <= x_pts[i + 1]:
                        t = (x_val - x_pts[i]) / (x_pts[i + 1] - x_pts[i])
                        return y_vals[i] + t * (y_vals[i + 1] - y_vals[i])
                return y_vals[-1] if x_val >= x_pts[-1] else y_vals[0]
            
            current_approx = new_approx
    
    return results

x0, y0 = 1.0, 1.0
a, b = 1.0, 2.0
x_points = [1.0, 1.2, 1.4, 1.6, 1.8, 2.0]

results = solve_picard(x0, y0, x_points, max_iterations=4)

exact_solution = [x**2 for x in x_points]
print("x\tТочное y(x)\tПриближение 4")
for i, x in enumerate(x_points):
    print(f"{x:.1f}\t{exact_solution[i]:.6f}\t{results[4][i]:.6f}")