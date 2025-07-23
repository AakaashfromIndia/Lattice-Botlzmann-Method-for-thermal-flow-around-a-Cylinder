import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.widgets import Slider
try:
    from numba import njit
    USE_NUMBA = True
except ImportError:
    USE_NUMBA = False

NX, NY = 420, 180
OB_RADIUS = NY // 9
OB_X, OB_Y = NX // 4, NY // 2

u_in = 0.08
omega = 1.0

c = np.array([
    [0, 0], [1, 0], [0, 1], [-1, 0], [0, -1],
    [1, 1], [-1, 1], [-1, -1], [1, -1]
])
w = np.array([4/9] + [1/9]*4 + [1/36]*4)

def initialize_fields():
    f = np.ones((NX, NY, 9))
    X, Y = np.ogrid[0:NX, 0:NY]
    obstacle = ((X - OB_X)**2 + (Y - OB_Y)**2) < OB_RADIUS**2
    return f, obstacle

def equilibrium_py(rho, ux, uy):
    feq = np.zeros((NX, NY, 9))
    u_sqr = ux**2 + uy**2
    for k in range(9):
        cu = 3 * (c[k, 0] * ux + c[k, 1] * uy)
        feq[:, :, k] = w[k] * rho * (1 + cu + 0.5 * cu**2 - 1.5 * u_sqr)
    return feq

def compute_macros_py(f, obstacle):
    rho = np.sum(f, axis=2)
    ux = np.sum(f * c[:, 0], axis=2) / rho
    uy = np.sum(f * c[:, 1], axis=2) / rho
    ux[obstacle] = 0
    uy[obstacle] = 0
    return rho, ux, uy

def bounce_back_py(f, obstacle):
    bb_pairs = [0, 3, 4, 1, 2, 7, 8, 5, 6]
    f_bb = f.copy()
    for k in range(9):
        f_bb[obstacle, k] = f[obstacle, bb_pairs[k]]
    return f_bb

def stream_py(f):
    f_streamed = np.empty_like(f)
    for k in range(9):
        f_streamed[:, :, k] = np.roll(
            np.roll(f[:, :, k], c[k, 0], axis=0),
                                c[k, 1], axis=1)
    return f_streamed

if USE_NUMBA:
    equilibrium = njit(equilibrium_py)
    compute_macros = njit(compute_macros_py)
    bounce_back = njit(bounce_back_py)
    stream = njit(stream_py)
else:
    equilibrium = equilibrium_py
    compute_macros = compute_macros_py
    bounce_back = bounce_back_py
    stream = stream_py

def inlet_bc_zouhe(f, rho, u):
    ux = u
    uy = 0.0
    rho_in = (f[0, :, 0] + f[0, :, 2] + f[0, :, 4] + 2*(f[0, :, 3] + f[0, :, 6] + f[0, :, 7])) / (1 - ux)
    f[0, :, 1] = f[0, :, 3] + (2/3)*rho_in*ux
    f[0, :, 5] = f[0, :, 7] + 0.5*(f[0, :, 4] - f[0, :, 2]) + (1/6)*rho_in*ux
    f[0, :, 8] = f[0, :, 6] + 0.5*(f[0, :, 2] - f[0, :, 4]) + (1/6)*rho_in*ux
    return f

def outlet_bc_zouhe(f):
    f[-1, :, 3] = f[-2, :, 3]
    f[-1, :, 6] = f[-2, :, 6]
    f[-1, :, 7] = f[-2, :, 7]
    return f

f, obstacle = initialize_fields()
rho, ux, uy = compute_macros(f, obstacle)

fig, ax = plt.subplots(figsize=(10, 4.5))
plt.subplots_adjust(left=0.25, bottom=0.25)
velocity_img = ax.imshow(np.sqrt(ux.T**2 + uy.T**2), origin='lower', cmap='jet',
                         vmin=0, vmax=u_in*2)
ax.set_title("Lattice Botlzmann Method for thermal flow around a Cylinder - Velocity Magnitude")

obstacle_img = ax.imshow(obstacle.T, origin='lower', cmap='binary', alpha=0.22, vmin=0, vmax=1)
cbar = plt.colorbar(velocity_img, ax=ax, fraction=0.032, pad=0.2)
cbar.set_label('Velocity Magnitude')

ax_u = plt.axes([0.25, 0.11, 0.6, 0.03])
slider_u = Slider(ax_u, 'Inlet Velocity', 0.01, 0.18, valinit=u_in, valstep=0.005)
ax_omega = plt.axes([0.25, 0.07, 0.6, 0.03])
slider_omega = Slider(ax_omega, 'Omega', 0.7, 1.99, valinit=omega, valstep=0.01)

def update(frame):
    global f, obstacle, rho, ux, uy, omega, u_in
    omega = slider_omega.val
    u_in = slider_u.val
    for substep in range(2):
        rho, ux, uy = compute_macros(f, obstacle)
        feq = equilibrium(rho, ux, uy)
        f += -omega * (f - feq)
        f = bounce_back(f, obstacle)
        f = stream(f)
        f = inlet_bc_zouhe(f, rho, u_in)
        f = outlet_bc_zouhe(f)
    rho, ux, uy = compute_macros(f, obstacle)
    velocity = np.sqrt(ux ** 2 + uy ** 2)
    velocity_img.set_array(velocity.T)
    return [velocity_img, obstacle_img]

ani = FuncAnimation(fig, update, frames=2000, interval=20, blit=False)
plt.show()
