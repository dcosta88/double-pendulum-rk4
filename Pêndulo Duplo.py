import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np

#Constantes
g = 9.81 #m s^(-2)
dt = 1/60
passos_por_frame = 4
dt_fisica = dt/passos_por_frame
#Condições do sistema
l1, l2 = 1, 1 #m
L = l1 + l2
m1, m2 = 1, 1 #kg

Y = np.array([np.pi, np.pi/2, 0, 0]) #[theta1, theta2, omega1, omega2]

historico_x2 = []
historico_y2 = []

#Funções
def energia(estado):
    theta1, theta2, omega1, omega2 = estado

    T1 = (m1 * l1 ** 2 * omega1 ** 2) / 2
    U1 = m1 * g * (l1 * (1 - np.cos(theta1)) + l2)
        
    T2 = (m2 * (l1 ** 2 * omega1 ** 2 + l2 ** 2 * omega2 ** 2 + 2 * l1* l2 * omega1 * omega2 * np.cos(theta1 - theta2))) / 2
    U2 = m2 * g * ((l1 * (1 - np.cos(theta1)) + l2) - l2 * np.cos(theta2))

    return T1 + U1 + T2 + U2

E0 = energia(Y)

def derivadas(Y):
    theta1, theta2, omega1, omega2 = Y
    d_theta = theta1 - theta2
    
    #Matriz M
    M = np.array([
            [(m1 + m2) * l1,          m2 * l2 * np.cos(d_theta)],
            [l1 * np.cos(d_theta),     l2]
        ])
    
    #Matriz F
    F = np.array([-g * (m1+m2) * np.sin(theta1) - m2 * l2 * np.sin(d_theta) * omega2 ** 2, 
                      l1 * np.sin(d_theta) * omega1 ** 2 - g * np.sin(theta2)])
    
    #Resolução do sistema
    alpha1, alpha2 = np.linalg.solve(M, F)

    return np.array([omega1, omega2, alpha1, alpha2])

def update(frame):
    global Y

    theta1, theta2, omega1, omega2 = Y

    for _ in range(passos_por_frame):
        k1 = derivadas(Y)
        k2 = derivadas(Y + k1 * dt_fisica / 2)
        k3 = derivadas(Y + k2 * dt_fisica / 2)
        k4 = derivadas(Y + k3 * dt_fisica)
        Y += (dt_fisica / 6) * (k1 + 2 * k2 + 2 * k3 + k4)

    #Posições e velocidades angulares
    theta1, theta2, omega1, omega2 = Y

    #Posições
    x1 = l1 * np.sin(theta1)
    y1 = l1 * (1 - np.cos(theta1)) + l2

    x2 = x1 + l2 * np.sin(theta2)
    y2 = y1 - l2 * np.cos(theta2)

    #Energia  
    E = energia(Y)
    delta_E = E - E0

    #Rasto
    historico_x2.append(x2)
    historico_y2.append(y2)

    if len(historico_x2) > 150:
        historico_x2.pop(0)
        historico_y2.pop(0)

    #Atualizações
    fio1.set_data([0, x1], [L, y1])
    fio2.set_data([x1, x2], [y1, y2])
    rasto.set_data(historico_x2, historico_y2)
    texto_dados.set_text(f"Energia Total: {E:.2f} J\nΔE: {delta_E:.5f} J")

    return fio1, fio2, rasto, texto_dados

#Animação
fig, ax = plt.subplots(figsize=(7, 7))
ax.set_xlim(-1.5 * L, 1.5 * L)
ax.set_ylim(-0.5 * L, 2.5 * L)
ax.set_aspect("equal")
fio1, = ax.plot([], [], "o-", color="black", markerfacecolor="red", markeredgecolor="black")
fio2, = ax.plot([], [], "o-", color="black", markerfacecolor="red", markeredgecolor="black")
rasto, = ax.plot([], [], "-", color="purple", alpha=0.5, linewidth=1)
texto_dados = ax.text(L * -1.4, L * 2.2, "", fontsize=10, fontfamily="monospace", 
                      bbox=dict(facecolor="pink", alpha=0.7, edgecolor="gray"))

ani = FuncAnimation(fig, update, interval=(dt * 10 ** 3), blit=True)

plt.xlabel("x (m)")
plt.ylabel("y (m)")
plt.title("Trajetória de um pêndulo duplo")
plt.grid(True)
plt.show()