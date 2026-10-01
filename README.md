# Double Pendulum Simulation

A Python simulation modeling the nonlinear chaotic motion of a double pendulum using the Lagrangian formalism, integrated via a custom 4th-order Runge-Kutta (RK4) algorithm with temporal sub-stepping.

![Simulation Demo](animation.mp4)

---

## Theoretical Background

The system consists of two coupled pendulums with point masses $m_1, m_2$ suspended by rigid, massless rods of lengths $l_1, l_2$. 

The equations of motion are derived using the **Euler-Lagrange equations**:

$$\frac{d}{dt}\left(\frac{\partial L}{\partial \dot{\theta}_i}\right) - \frac{\partial L}{\partial \theta_i} = 0 \quad (i = 1, 2)$$

where the Lagrangian is $L = T - U$. Due to the nonlinear coupling, the equations are arranged into a matrix system:

$$M(\theta) \begin{bmatrix} \ddot{\theta}_1 \\ \ddot{\theta}_2 \end{bmatrix} = F(\theta, \dot{\theta})$$

At each numerical step, angular accelerations $(\ddot{\theta}_1, \ddot{\theta}_2)$ are obtained by solving this linear system.

---

## Computational Implementation

- **Explicit RK4 Integrator:** Solves the first-order system of four coupled differential equations. Unlike lower-order methods (e.g., explicit Euler), RK4 achieves local truncation error of $\mathcal{O}(dt^5)$ and global error of $\mathcal{O}(dt^4)$, preserving physical trajectory stability.
- **Sub-stepping:** Runs multiple numerical integration steps per animation frame ($dt_{\text{physics}} = dt / N$) to prevent numerical divergence during high angular velocity regimes.
- **Energy Conservation Tracking ($\Delta E$):** Computes total mechanical energy ($E = T + U$) in real time. The drift $\Delta E = E(t) - E_0$ serves as a direct quantitative benchmark of solver fidelity.
- **Visualization:** Built with `matplotlib.animation.FuncAnimation` utilizing `blit=True` for high-performance rendering, persistent trajectory tracing, and real-time metrics display.

---

## Requirements & Installation

Requires Python 3.8+ with standard scientific libraries:

```bash
pip install numpy matplotlib
