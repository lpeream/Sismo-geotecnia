# -*- coding: utf-8 -*-
"""
Paso 3: genera los 3 acelerogramas sinteticos finales, compatibles con
el espectro NSR-10 de Neiva (perfil B, roca), segun A.2.7 del Reglamento.
"""

from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

from generar_acelerogramas import espectro_respuesta
from ajuste_espectral import sa_nsr10, generar_semilla, ajustar_espectro, T0, TC, TL

ROOT = Path(__file__).resolve().parents[1]
FIGURAS = ROOT / "figuras"
DATOS = ROOT / "datos" / "datos depurados"

dt = 0.01
duracion = 20.0
periodos_ajuste = np.geomspace(0.05, 3.5, 25)
sa_objetivo = sa_nsr10(periodos_ajuste)
periodos_verificacion = np.geomspace(0.02, 4.0, 60)
sa_target_verif = sa_nsr10(periodos_verificacion)

acelerogramas = []
espectros_finales = []

for numero, semilla in enumerate([1, 7, 42], start=1):
    print(f"\nGenerando acelerograma {numero} (semilla={semilla})...")
    t, acel_inicial = generar_semilla(dt, duracion, semilla_rng=semilla)
    acel_final = ajustar_espectro(acel_inicial, dt, periodos_ajuste, sa_objetivo, n_iter=12)
    acelerogramas.append(acel_final)

    sa_final = espectro_respuesta(acel_final, dt, periodos_verificacion)
    espectros_finales.append(sa_final)

    # Guardar el acelerograma como archivo de texto (tiempo, aceleracion en g)
    salida = DATOS.parent / "acelerogramas_sinteticos"
    salida.mkdir(exist_ok=True)
    archivo = salida / f"acelerograma_neiva_roca_{numero}.txt"
    np.savetxt(archivo, np.column_stack([t, acel_final]),
               header="tiempo_s  aceleracion_g", fmt="%.4f")
    print(f"  -> Guardado: {archivo}")

espectro_promedio = np.mean(espectros_finales, axis=0)

# --- Verificacion final: el promedio de los 3 no debe caer por debajo del
# objetivo en el rango 0.2T a 1.5T (criterio A.2.7(c), aplicado de forma
# general a todo el rango de periodos relevantes para este estudio) ---
razon_promedio = espectro_promedio / sa_target_verif
print(f"\nRazon promedio/objetivo: min={razon_promedio.min():.2f}  "
      f"max={razon_promedio.max():.2f}  media={razon_promedio.mean():.2f}")

# --- Figura final ---
fig, axes = plt.subplots(2, 1, figsize=(11, 8))
colores = ["#1f5d75", "#7a3b3b", "#3b7a4e"]
t_plot = np.arange(0, duracion, dt)
for i, acel in enumerate(acelerogramas):
    axes[0].plot(t_plot, acel + i * 1.0, color=colores[i], linewidth=0.5,
                 label=f"Acelerograma {i+1}")
axes[0].set_xlabel("Tiempo (s)")
axes[0].set_ylabel("Aceleracion (g) [desplazados para claridad]")
axes[0].set_title("3 acelerogramas sinteticos compatibles - Neiva, perfil B (roca)")
axes[0].legend(loc="upper right", fontsize=8)
axes[0].grid(alpha=0.3)

axes[1].plot(periodos_verificacion, sa_target_verif, color="black", linewidth=2.5, label="NSR-10 objetivo")
for i, sa in enumerate(espectros_finales):
    axes[1].plot(periodos_verificacion, sa, color=colores[i], linewidth=1, alpha=0.6,
                 linestyle="--", label=f"Acelerograma {i+1}")
axes[1].plot(periodos_verificacion, espectro_promedio, color="#b24732", linewidth=2,
             label="Promedio de los 3")
axes[1].set_xlabel("Periodo T (s)")
axes[1].set_ylabel("Sa (g)")
axes[1].set_title("Espectros de respuesta individuales, promedio, y objetivo NSR-10")
axes[1].legend(fontsize=8)
axes[1].grid(alpha=0.3)

fig.tight_layout()
fig.savefig(FIGURAS / "acelerogramas_finales_neiva.png", dpi=220, bbox_inches="tight")
print("\n-> Figura guardada en figuras/acelerogramas_finales_neiva.png")