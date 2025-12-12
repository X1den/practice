import math

a, b, n = 1.0, 2.0, 12

alpha, beta, gamma, delta = 1.0, -5.0, 2.0, 1.0

A, B = -1.0, 4.0

def p(x):
    return 2.0 / math.sqrt(x)

def q(x):
    return -x

def f(x):
    return 1.0 / x

def F(x, y):
    u, up = y
    upp = f(x) - p(x)*up - q(x)*u
    return [up, upp]

def rk4(y0, x0, x1, nsteps):
    h = (x1 - x0) / nsteps
    xs = [x0]
    ys = [y0[:]]

    for i in range(nsteps):
        x = xs[-1]
        y = ys[-1]

        k1 = F(x, y)
        k2 = F(x + h/2, [y[j] + h*k1[j]/2 for j in range(2)])
        k3 = F(x + h/2, [y[j] + h*k2[j]/2 for j in range(2)])
        k4 = F(x + h,   [y[j] + h*k3[j]   for j in range(2)])

        y_new = [
            y[j] + (h/6)*(k1[j] + 2*k2[j] + 2*k3[j] + k4[j])
            for j in range(2)
        ]

        xs.append(x + h)
        ys.append(y_new)

    return xs, ys

v0 = [0.0, A/alpha]

w0 = [1.0, -beta/alpha]

xs, vs = rk4(v0, a, b, n)
_, ws = rk4(w0, a, b, n)

v  = [row[0] for row in vs]
vp = [row[1] for row in vs]

w  = [row[0] for row in ws]
wp = [row[1] for row in ws]

numer = B - delta*vp[-1] - gamma*v[-1]
denom = delta*wp[-1] + gamma*w[-1]
c = numer / denom

u  = [v[i]  + c*w[i]  for i in range(len(v))]
up = [vp[i] + c*wp[i] for i in range(len(v))]

print("c =", c)
print("x          u(x)           u'(x)")
for i in range(len(xs)):
    print(f"{xs[i]:.6f}   {u[i]:.8f}   {up[i]:.8f}")