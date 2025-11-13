import math

def f(x): 
    return x**3 * math.exp(-0.12 * x**1.5)

def cubic_spline_M(x_nodes, y_nodes):
    n = len(x_nodes) - 1
    h = [x_nodes[i+1]-x_nodes[i] for i in range(n)]
    A = [[0]*(n+1) for _ in range(n+1)]
    B = [0]*(n+1)
    A[0][0] = A[n][n] = 1
    for i in range(1,n):
        A[i][i-1] = h[i-1]
        A[i][i] = 2*(h[i-1]+h[i])
        A[i][i+1] = h[i]
        B[i] = 6*((y_nodes[i+1]-y_nodes[i])/h[i] - (y_nodes[i]-y_nodes[i-1])/h[i-1])
    
    for i in range(n+1):
        factor = A[i][i]
        for j in range(i, n+1):
            A[i][j] /= factor
        B[i] /= factor
        for k in range(i+1, n+1):
            factor = A[k][i]
            for j in range(i, n+1):
                A[k][j] -= factor*A[i][j]
            B[k] -= factor*B[i]
    M = [0]*(n+1)
    for i in range(n, -1, -1):
        M[i] = B[i]
        for j in range(i+1, n+1):
            M[i] -= A[i][j]*M[j]
    return M

def linspace(a, b, n):
    h = (b-a)/(n-1)
    return [a + i*h for i in range(n)]

def spline_integral(f, a, b, n):
    x_nodes = linspace(a,b,n)
    y_nodes = [f(x) for x in x_nodes]
    M = cubic_spline_M(x_nodes, y_nodes)
    h = [x_nodes[i+1]-x_nodes[i] for i in range(n-1)]
    integral = 0
    for i in range(n-1):
        integral += (1/2)*h[i]*(y_nodes[i]+y_nodes[i+1]) - (1/24)*(h[i]**3)*(M[i]+M[i+1])
    return integral, x_nodes, y_nodes, M

def spline_error(f, a, b, n):
    I_spline, x_nodes, y_nodes, M = spline_integral(f, a, b, n)
    x_dense = linspace(a,b,1000)
    S_vals = []
    for x in x_dense:
        for i in range(len(x_nodes)-1):
            if x_nodes[i]<=x<=x_nodes[i+1]:
                h_i = x_nodes[i+1]-x_nodes[i]
                term = (M[i]*(x_nodes[i+1]-x)**3 + M[i+1]*(x-x_nodes[i])**3)/(6*h_i)
                term += (y_nodes[i]-M[i]*h_i**2/6)*(x_nodes[i+1]-x)/h_i
                term += (y_nodes[i+1]-M[i+1]*h_i**2/6)*(x-x_nodes[i])/h_i
                S_vals.append(term)
                break
    max_error = max(abs(f(x)-Sx) for x,Sx in zip(x_dense,S_vals))
    return I_spline, max_error

a, b = 1.5, 3.75
n = 16
I_spline, max_err = spline_error(f, a, b, n)
print(f"Интеграл сплайном: {I_spline:.10f}, макс. погрешность |f-S| ≈ {max_err:.2e}")