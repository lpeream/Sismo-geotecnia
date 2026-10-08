# -*- coding: utf-8 -*-
"""
Paso 1: funcion de espectro de respuesta (metodo de Newmark-beta, forma
incremental estandar), con una prueba de verificacion antes de usarla en
algo mas complejo.
"""

import numpy as np


def espectro_respuesta(acelerograma, dt, periodos, amortiguamiento=0.05):
    """
    Calcula el espectro de pseudo-aceleracion Sa(T) para un oscilador de
    1 grado de libertad, masa unitaria, usando Newmark-beta (promedio
    constante: beta=1/4, gamma=1/2), forma incremental estandar.
    """
    beta, gamma = 0.25, 0.5
    n = len(acelerograma)
    Sa = np.zeros(len(periodos))
    m = 1.0

    for idx, T in enumerate(periodos):
        omega = 2 * np.pi / T
        k = omega ** 2 * m
        c = 2 * amortiguamiento * omega * m

        u = np.zeros(n)
        v = np.zeros(n)
        a = np.zeros(n)
        a[0] = (-m * acelerograma[0] - c * v[0] - k * u[0]) / m

        k_hat = k + (gamma / (beta * dt)) * c + (1.0 / (beta * dt ** 2)) * m
        a1 = (1.0 / (beta * dt)) * m + (gamma / beta) * c
        a2 = (1.0 / (2 * beta)) * m + dt * (gamma / (2 * beta) - 1) * c

        for i in range(n - 1):
            delta_p = -m * (acelerograma[i + 1] - acelerograma[i])
            delta_p_hat = delta_p + a1 * v[i] + a2 * a[i]
            delta_u = delta_p_hat / k_hat
            delta_v = (gamma / (beta * dt)) * delta_u - (gamma / beta) * v[i] + dt * (1 - gamma / (2 * beta)) * a[i]
            delta_a = (1.0 / (beta * dt ** 2)) * delta_u - (1.0 / (beta * dt)) * v[i] - (1.0 / (2 * beta)) * a[i]

            u[i + 1] = u[i] + delta_u
            v[i + 1] = v[i] + delta_v
            a[i + 1] = a[i] + delta_a

        Sa[idx] = np.max(np.abs(omega ** 2 * u))  # pseudo-aceleracion

    return Sa


if __name__ == "__main__":
    dt = 0.005
    t = np.arange(0, 5, dt)
    frecuencia_prueba = 2.0  # Hz -> periodo de 0.5 s
    acel_prueba = np.zeros_like(t)
    duracion_pulso = 1.0 / frecuencia_prueba
    mascara = t < duracion_pulso
    acel_prueba[mascara] = np.sin(2 * np.pi * frecuencia_prueba * t[mascara])

    periodos_prueba = np.array([0.1, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0])
    Sa_prueba = espectro_respuesta(acel_prueba, dt, periodos_prueba)

    print("Prueba de verificacion: pulso a 0.5s de periodo")
    for T, sa in zip(periodos_prueba, Sa_prueba):
        print(f"  T={T:.2f}s  Sa={sa:.4f}")
    print("\nSe espera que el valor mas alto este cerca de T=0.5s (el periodo del pulso).")