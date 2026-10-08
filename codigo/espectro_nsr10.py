# -*- coding: utf-8 -*-
"""
Espectro elastico de aceleraciones de diseno, NSR-10, Capitulo A.2.
Ciudad: Neiva (Huila). Perfil de suelo: B (roca de rigidez media).
Valores oficiales tomados de la Tabla A.2.3-2 de la NSR-10 (Titulo A).
"""

from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
FIGURAS = ROOT / "figuras"
DATOS = ROOT / "datos" / "datos depurados"
FIGURAS.mkdir(exist_ok=True)

# --- Parametros de entrada (NSR-10, Tabla A.2.3-2 y A.2.4-3/A.2.4-4) ---
Aa = 0.25   # Neiva, aceleracion horizontal pico efectiva
Av = 0.25   # Neiva, velocidad horizontal pico efectiva
Fa = 1.0    # Perfil B (roca), independiente de la intensidad
Fv = 1.0    # Perfil B (roca), independiente de la intensidad
I = 1.0     # Coeficiente de importancia, grupo de uso I (generico)

# --- Periodos caracteristicos (Ec. A.2.6-2 y A.2.6-4) ---
TC = 0.48 * (Av * Fv) / (Aa * Fa)
TL = 2.4 * Fv
T0 = 0.1 * (Av * Fv) / (Aa * Fa)

print(f"Aa={Aa}  Av={Av}  Fa={Fa}  Fv={Fv}  I={I}")
print(f"T0 = {T0:.3f} s")
print(f"TC = {TC:.3f} s")
print(f"TL = {TL:.3f} s")


def sa_nsr10(T):
    """Espectro elastico de aceleraciones Sa(T), fraccion de g. Ec. A.2.6-1 a A.2.6-5."""
    T = np.atleast_1d(T).astype(float)
    Sa = np.zeros_like(T)

    # Zona 1: T < T0 (Ec. A.2.6-7, usada tambien como transicion inicial)
    mascara_1 = T < T0
    Sa[mascara_1] = 2.5 * Aa * Fa * I * (0.4 + 0.6 * T[mascara_1] / T0)

    # Zona 2: T0 <= T <= TC (aceleracion constante, Ec. A.2.6-3)
    mascara_2 = (T >= T0) & (T <= TC)
    Sa[mascara_2] = 2.5 * Aa * Fa * I

    # Zona 3: TC < T <= TL (Ec. A.2.6-1)
    mascara_3 = (T > TC) & (T <= TL)
    Sa[mascara_3] = 1.2 * Av * Fv * I / T[mascara_3]

    # Zona 4: T > TL (Ec. A.2.6-5)
    mascara_4 = T > TL
    Sa[mascara_4] = 1.2 * Av * Fv * TL * I / T[mascara_4] ** 2

    return Sa


# --- Genera el espectro completo ---
periodos = np.linspace(0.01, 4.0, 800)
Sa = sa_nsr10(periodos)

# --- Grafica ---
fig, ax = plt.subplots(figsize=(9, 5.5))
ax.plot(periodos, Sa, color="#1f5d75", linewidth=2)
ax.axvline(T0, color="gray", linestyle=":", linewidth=1)
ax.axvline(TC, color="gray", linestyle=":", linewidth=1)
ax.axvline(TL, color="gray", linestyle=":", linewidth=1)
ax.text(T0, ax.get_ylim()[1]*0.02, "$T_0$", ha="center", fontsize=9)
ax.text(TC, ax.get_ylim()[1]*0.02, "$T_C$", ha="center", fontsize=9)
ax.text(TL, ax.get_ylim()[1]*0.02, "$T_L$", ha="center", fontsize=9)
ax.set_xlabel("Periodo T (s)")
ax.set_ylabel(r"$S_a$ (fraccion de g)")
ax.set_title(f"Espectro elastico de diseno NSR-10 - Neiva, perfil B (roca)\n"
             f"Aa={Aa}, Av={Av}, Fa={Fa}, Fv={Fv}, I={I}")
ax.grid(alpha=0.3)
fig.tight_layout()
fig.savefig(FIGURAS / "espectro_nsr10_neiva_roca.png", dpi=220, bbox_inches="tight")
print("\n-> Figura guardada en figuras/espectro_nsr10_neiva_roca.png")

# --- Tabla de puntos clave para el informe ---
puntos_clave = pd.DataFrame({
    "T (s)": [0.0, T0, TC, 1.0, TL, 3.0, 4.0],
})
puntos_clave["Sa (g)"] = sa_nsr10(puntos_clave["T (s)"].to_numpy())
print("\nPuntos clave del espectro:")
print(puntos_clave.to_string(index=False))
puntos_clave.to_csv(DATOS / "espectro_nsr10_puntos_clave.csv", index=False, encoding="utf-8-sig")