# -*- coding: utf-8 -*-
"""
Paso 2: genera un acelerograma sintetico y ajusta su contenido de
frecuencias iterativamente para que su espectro de respuesta converja
al espectro NSR-10 objetivo (Neiva, perfil B).
"""

from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

from generar_acelerogramas import espectro_respuesta

ROOT = Path(__file__).resolve().parents[1]
FIGURAS = ROOT / "figuras"

# --- Espectro objetivo NSR-10 (mismos parametros que espectro_nsr10.py) ---
Aa, Av, Fa, Fv, I = 0.25, 0.25, 1.0, 1.0, 1.0
T0 = 0.1 * (Av * Fv) / (Aa * Fa)
TC = 0.48 * (Av * Fv) / (Aa * Fa)
TL = 2.4 * Fv


def sa_nsr10(T):
    T = np.atleast_1d(T).astype(float)
    Sa = np.zeros_like(T)
    m1 = T < T0
    Sa[m1] = 2.5 * Aa * Fa * I * (0.4 + 0.6 * T[m1] / T0)
    m2 = (T >= T0) & (T <= TC)
    Sa[m2] = 2.5 * Aa * Fa * I
    m3 = (T > TC) & (T <= TL)
    Sa[m3] = 1.2 * Av * Fv * I / T[m3]
    m4 = T > TL
    Sa[m4] = 1.2 * Av * Fv * TL * I / T[m4] ** 2
    return Sa


def envolvente_sismo(t, t_subida=2.0, t_fuerte=8.0, t_total=20.0):
    """Envolvente tipo Jennings: sube, se mantiene fuerte, decae exponencial."""
    env = np.ones_like(t)
    subida = t < t_subida
    env[subida] = (t[subida] / t_subida) ** 2
    decae = t > (t_subida + t_fuerte)
    td = t[decae] - (t_subida + t_fuerte)
    env[decae] = np.exp(-0.15 * td)
    return env


def generar_semilla(dt, duracion, semilla_rng):
    t = np.arange(0, duracion, dt)
    rng = np.random.default_rng(semilla_rng)
    ruido = rng.standard_normal(len(t))
    env = envolvente_sismo(t)
    return t, ruido * env


def ajustar_espectro(acel, dt, periodos_ajuste, sa_objetivo, n_iter=12):
    """Ajuste iterativo en el dominio de la frecuencia (metodo simplificado)."""
    n = len(acel)
    for iteracion in range(n_iter):
        sa_actual = espectro_respuesta(acel, dt, periodos_ajuste)
        razon = sa_objetivo / np.maximum(sa_actual, 1e-6)
        razon = np.clip(razon, 0.5, 2.0)  # limitar cambios bruscos por iteracion

        freqs_ajuste = 1.0 / periodos_ajuste  # Hz
        orden = np.argsort(freqs_ajuste)

        Y = np.fft.rfft(acel)
        freqs_fft = np.fft.rfftfreq(n, dt)
        escala = np.interp(freqs_fft, freqs_ajuste[orden], razon[orden],
                            left=razon[orden][0], right=razon[orden][-1])
        Y_ajustado = Y * escala
        acel = np.fft.irfft(Y_ajustado, n=n)

        error = np.mean(np.abs(sa_objetivo - sa_actual) / sa_objetivo)
        print(f"  Iteracion {iteracion+1}: error medio relativo = {error*100:.1f}%")

    return acel


# --- Prueba con 1 acelerograma ---
dt = 0.01
duracion = 20.0
periodos_ajuste = np.geomspace(0.05, 3.5, 25)
sa_objetivo = sa_nsr10(periodos_ajuste)

t, acel_inicial = generar_semilla(dt, duracion, semilla_rng=1)
print("Ajustando acelerograma de prueba...")
acel_final = ajustar_espectro(acel_inicial, dt, periodos_ajuste, sa_objetivo)

# --- Verificacion final ---
periodos_verificacion = np.geomspace(0.02, 4.0, 60)
sa_final = espectro_respuesta(acel_final, dt, periodos_verificacion)
sa_target_verif = sa_nsr10(periodos_verificacion)

fig, axes = plt.subplots(1, 2, figsize=(13, 5))
axes[0].plot(t, acel_final, color="#1f5d75", linewidth=0.6)
axes[0].set_xlabel("Tiempo (s)")
axes[0].set_ylabel("Aceleracion (g)")
axes[0].set_title("Acelerograma sintetico de prueba")
axes[0].grid(alpha=0.3)

axes[1].plot(periodos_verificacion, sa_target_verif, color="black", linewidth=2, label="NSR-10 objetivo")
axes[1].plot(periodos_verificacion, sa_final, color="#b24732", linewidth=1.5, linestyle="--", label="Acelerograma ajustado")
axes[1].set_xlabel("Periodo T (s)")
axes[1].set_ylabel("Sa (g)")
axes[1].set_title("Verificacion del ajuste espectral")
axes[1].legend()
axes[1].grid(alpha=0.3)

fig.tight_layout()
fig.savefig(FIGURAS / "prueba_ajuste_espectral.png", dpi=200, bbox_inches="tight")
print("\n-> Figura guardada en figuras/prueba_ajuste_espectral.png")