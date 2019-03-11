import numpy as numpy
from matplotlib import pyplot, cm
from matplotlib.colors import Normalize
import matplotlib
matplotlib.use('Agg')


nx = 41
ny = 41
dx = 2/float(nx-1)
dy = 2/float(ny-1)

u = numpy.ones((ny, nx))
v = numpy.ones((ny, nx))

sigma = 0.001
nu = 0.01
dt = sigma*dx*dy/nu


def equation_of_motion(u, v, params):
    [dt, dx, dy, nu, nozzle_u, nozzle_v, nx, ny] = params
    # generate the next state as a func of old state
    # un = u.copy()
    # vn = v.copy()
    un = u[1:-1, 1:-1]
    ufc = u[1:-1, 2:]
    ubc = u[1:-1, 0:-2]
    ufr = u[2:, 1:-1]
    ubr = u[0:-2, 1:-1]
    # Setting y-velocity grid components
    vn = v[1:-1, 1:-1]
    vfc = v[1:-1, 2:]
    vbc = v[1:-1, 0:-2]
    vfr = v[2:, 1:-1]
    vbr = v[0:-2, 1:-1]
    # Solving for u_next
    u_num1 = ufr-(2*un)+ubr
    u_num2 = ufc-(2*un)+ubc
    u_first = nu*((u_num1/(dx**2)) + (u_num2/(dy**2)))
    u_second = -un*((un-ubr)/dx)
    u_third = -vn*((un-ubc)/dy)
    u_next = un + dt*(u_first+u_second+u_third)
    u[1:-1, 1:-1] = u_next
    # Solving for v_next
    v_num1 = vfr-(2*vn)+vbr
    v_num2 = vfc-(2*vn)+vbc
    v_first = nu*((v_num1/(dx**2)) + (v_num2/(dy**2)))
    v_second = -un*((vn-vbr)/dx)
    v_third = -vn*((vn-vbc)/dy)
    v_next = vn + dt*(v_first+v_second+v_third)
    v[1:-1, 1:-1] = v_next
    return (u, v)


def boundary(u, v, params, t_step):
    [dt, dx, dy, nu, nozzle_u, nozzle_v, nx, ny] = params
    print u
    u[0, :] = 0
    u[-1, :] = 0
    u[:, 0] = 0
    u[:, -1] = 0

    v[0, :] = 0
    v[-1, :] = 0
    v[:, 0] = 0
    v[:, -1] = 0
    # Special nozzle BC
    u[ny//2-2:ny//2+2, 0] = nozzle_u[t_step]
    v[ny//2-2:ny//2+2, 0] = nozzle_v[t_step]
    return (u, v)


def evolve(u, v, params, steps):
    for i in range(steps):
        (u, v) = equation_of_motion(u, v, params)
        (u, v) = boundary(u, v, params, i)
        if i % 50 == 0:
            ax = pyplot.figure()
            norm = Normalize()
            magnitude = numpy.sqrt(u[::2]**2 + v[::2]**2)
            pyplot.quiver(u[::2], v[::2], norm(magnitude), scale=50, cmap=pyplot.cm.jet)
            ax.savefig('frame'+str(i).zfill(5)+'.png', dpi=300)
            ax.clear()
    return (u, v)


nt = 2510  # number of steps simulating
# Assigning initial conditions
initial_u = numpy.zeros((nx, ny))
initial_v = numpy.zeros((nx, ny))
# Special BC for nozzle located at (0,1)
nozzle_u = numpy.append(10*numpy.ones(1000), numpy.zeros(nt))
nozzle_v = numpy.append(10*numpy.ones(1000), numpy.zeros(nt))

params = [dt, dx, dy, nu, nozzle_u, nozzle_v, nx, ny]

(final_u, final_v) = evolve(initial_u, initial_v, params, nt)
