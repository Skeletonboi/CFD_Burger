import numpy as np
import matplotlib.pyplot as plt

u = np.array(range(10))
un = u[1:-1]
ub = u[0:-2]
u_next = un - ub
print u_next


def step_u(u, params):
    [dt, dx, c] = params
    A = c*dt/dx
    B = dt
    u_next = np.zeros_like(u)
    un = u[1:-1]
    ub = u[0:-2]
    linear = (1+A)*un-A*ub
    nonlinear = -B*un**2
    u_next[1:-1] = linear + nonlinear
    return u_next


def BC_u(u, params):
    u[0] = 0
    u[-1] = 0
    return u


def sim(nt, u, params):
    for i in range(nt):
        u = step_u(u, params)
        u = BC_u(u, params)
    return u


nt = 100
L = 10
dt = 0.01
dx = 0.1
c = 1
params = [dt, dx, c]


def gauss(x):
    return np.exp(-x**2)


x = np.arange(0, L, dx)
u_init = np.array(list(map(gauss, x-np.mean(x))))

u_after = sim(nt, u_init, params)

plt.figure()
plt.plot(x, u_init)
plt.show()

plt.figure()
plt.plot(x, u_after)
plt.show()
