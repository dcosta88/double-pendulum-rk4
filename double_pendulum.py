import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np

# Constants
g = 9.81 # m s^(-2)
dt = 1/60
substeps = 4
dt_physics = dt / substeps

# System parameters
l1, l2 = 1, 1 # m
L = l1 + l2
m1, m2 = 1, 1 # kg

Y = np.array([np.pi, np.pi/2, 0, 0]) # [theta1, theta2, omega1, omega2]

history_x2 = []
history_y2 = []

# Functions
def energy(state):
    theta1, theta2, omega1, omega2 = state

    T1 = (m1 * l1 ** 2 * omega1 ** 2) / 2
    U1 = m1 * g * (l1 * (1 - np.cos(theta1)) + l2)
        
    T2 = (m2 * (l1 ** 2 * omega1 ** 2 + l2 ** 2 * omega2 ** 2 + 2 * l1 * l2 * omega1 * omega2 * np.cos(theta1 - theta2))) / 2
    U2 = m2 * g * ((l1 * (1 - np.cos(theta1)) + l2) - l2 * np.cos(theta2))

    return T1 + U1 + T2 + U2

E0 = energy(Y)

def derivatives(Y):
    theta1, theta2, omega1, omega2 = Y
    d_theta = theta1 - theta2
    
    # Matrix M
    M = np.array([
            [(m1 + m2) * l1,          m2 * l2 * np.cos(d_theta)],
            [l1 * np.cos(d_theta),     l2]
        ])
    
    # Vector F
    F = np.array([-g * (m1 + m2) * np.sin(theta1) - m2 * l2 * np.sin(d_theta) * omega2 ** 2, 
                  l1 * np.sin(d_theta) * omega1 ** 2 - g * np.sin(theta2)])
    
    # System resolution
    alpha1, alpha2 = np.linalg.solve(M, F)

    return np.array([omega1, omega2, alpha1, alpha2])

def update(frame):
    global Y

    theta1, theta2, omega1, omega2 = Y

    for _ in range(substeps):
        k1 = derivatives(Y)
        k2 = derivatives(Y + k1 * dt_physics / 2)
        k3 = derivatives(Y + k2 * dt_physics / 2)
        k4 = derivatives(Y + k3 * dt_physics)
        Y += (dt_physics / 6) * (k1 + 2 * k2 + 2 * k3 + k4)

    # Angular positions and velocities
    theta1, theta2, omega1, omega2 = Y

    # Cartesian coordinates
    x1 = l1 * np.sin(theta1)
    y1 = l1 * (1 - np.cos(theta1)) + l2

    x2 = x1 + l2 * np.sin(theta2)
    y2 = y1 - l2 * np.cos(theta2)

    # Energy calculations
    E = energy(Y)
    delta_E = E - E0

    # Trajectory trace
    history_x2.append(x2)
    history_y2.append(y2)

    if len(history_x2) > 150:
        history_x2.pop(0)
        history_y2.pop(0)

    # Visual updates
    rod1.set_data([0, x1], [L, y1])
    rod2.set_data([x1, x2], [y1, y2])
    trace.set_data(history_x2, history_y2)
    data_text.set_text(f"Total Energy: {E:.2f} J\nΔE: {delta_E:.5f} J")

    return rod1, rod2, trace, data_text

# Animation setup
fig, ax = plt.subplots(figsize=(7, 7))
ax.set_xlim(-1.5 * L, 1.5 * L)
ax.set_ylim(-0.5 * L, 2.5 * L)
ax.set_aspect("equal")
rod1, = ax.plot([], [], "o-", color="black", markerfacecolor="red", markeredgecolor="black")
rod2, = ax.plot([], [], "o-", color="black", markerfacecolor="red", markeredgecolor="black")
trace, = ax.plot([], [], "-", color="purple", alpha=0.5, linewidth=1)
data_text = ax.text(L * -1.4, L * 2.2, "", fontsize=10, fontfamily="monospace", 
                    bbox=dict(facecolor="pink", alpha=0.7, edgecolor="gray"))

ani = FuncAnimation(fig, update, interval=(dt * 10 ** 3), blit=True)

plt.xlabel("x (m)")
plt.ylabel("y (m)")
plt.title("Double Pendulum Trajectory")
plt.grid(True)
plt.show()
