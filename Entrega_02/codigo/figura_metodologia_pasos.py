"""Figura de la metodología: los cinco pasos del análisis con datos del estudio.

Paso 1 sismofuentes, Paso 2 recurrencia, Paso 3 atenuación, Paso 4 amenaza y
Paso 5 espectro de aceleración. Lee los resultados ya calculados por los
demás scripts; no recalcula nada.
"""
from __future__ import annotations

import importlib.util
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

DELIVERY_ROOT = Path(__file__).resolve().parents[1]
CATALOG_PATH = (DELIVERY_ROOT / "recursos" / "Qgis" / "datos"
                / "Sismos_260_clasificados_20261004_192107.csv")
STATS = DELIVERY_ROOT / "datos" / "datos estadistica"
HAZARD_PATH = STATS / "integracion paso4" / "curva_amenaza_neiva_roca_paso4.csv"
SPECTRUM_PATH = STATS / "espectros paso5" / "espectro_diseno_roca_paso5.csv"
PARAMETERS_PATH = STATS / "espectros paso5" / "parametros_espectro_diseno_paso5.csv"
SPECTRA_SCRIPT = DELIVERY_ROOT / "codigo" / "espectros_neiva.py"
FIGURE_DIRECTORY = DELIVERY_ROOT / "figuras"

SITE_LAT, SITE_LON = 2.9262767, -75.2892111
SOURCES = ("Cortical", "Subducción interplaca", "Subducción intraplaca")
SHORT = {"Cortical": "Cortical", "Subducción interplaca": "Interplaca",
         "Subducción intraplaca": "Intraplaca"}
COLORS = {"Cortical": "#1f77b4", "Subducción interplaca": "#ff7f0e",
          "Subducción intraplaca": "#2ca02c"}
# Parámetros Gutenberg-Richter (Tabla 2) y magnitudes máximas adoptadas.
GR = {"Cortical": (4.197, 0.882, 7.55), "Subducción interplaca": (3.911, 0.964, 6.60),
      "Subducción intraplaca": (4.025, 0.821, 7.90)}
# Profundidad típica de cada fuente para las curvas de atenuación ilustrativas.
DEPTH_KM = {"Cortical": 10.0, "Subducción interplaca": 20.0, "Subducción intraplaca": 95.0}


def _load_spectra_module():
    spec = importlib.util.spec_from_file_location("espectros", SPECTRA_SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main() -> None:
    spectra = _load_spectra_module()
    catalog = pd.read_csv(CATALOG_PATH, encoding="latin1")
    catalog["fuente"] = catalog["fuente"].str.replace("�", "ó")
    hazard = pd.read_csv(HAZARD_PATH, encoding="utf-8-sig")
    hazard.columns = [c.replace("�", "ó") for c in hazard.columns]
    design = pd.read_csv(SPECTRUM_PATH, encoding="utf-8-sig")
    params = pd.read_csv(PARAMETERS_PATH, encoding="utf-8-sig").iloc[0]

    plt.rcParams.update({"font.size": 10})
    fig, axes = plt.subplots(1, 5, figsize=(17, 3.9), dpi=300)
    titles = ("Paso 1\nSismofuentes", "Paso 2\nRecurrencia", "Paso 3\nAtenuación",
              "Paso 4\nAmenaza", "Paso 5\nEspectro de aceleración")
    for ax, title in zip(axes, titles):
        ax.set_title(title, fontsize=11, fontweight="bold")
        ax.grid(True, alpha=0.25, linewidth=0.6)

    # Paso 1: sismicidad clasificada alrededor del sitio.
    ax = axes[0]
    for source in SOURCES:
        sub = catalog[catalog["fuente"] == source]
        ax.scatter(sub["longitud"], sub["latitud"], s=3, color=COLORS[source],
                   alpha=0.6, linewidths=0, label=SHORT[source])
    ax.plot(SITE_LON, SITE_LAT, marker="*", color="black", markersize=11)
    ax.annotate("Neiva", (SITE_LON, SITE_LAT), xytext=(4, 4),
                textcoords="offset points", fontsize=9)
    ax.set_xlabel("Longitud (°)")
    ax.set_ylabel("Latitud (°)")
    ax.set_aspect("equal", adjustable="datalim")

    # Paso 2: rectas de Gutenberg-Richter truncadas en Mw_max.
    ax = axes[1]
    for source in SOURCES:
        a, b, m_max = GR[source]
        m = np.linspace(3.0, m_max, 100)
        ax.plot(m, 10 ** (a - b * m), color=COLORS[source], linewidth=2)
    ax.set_yscale("log")
    ax.set_xlabel("Magnitud Mw")
    ax.set_ylabel("λ (Mw ≥ M) (1/año)")

    # Paso 3: PGA mediana de Zhao et al. (2006) para Mw 6.5.
    ax = axes[2]
    for source in SOURCES:
        h = DEPTH_KM[source]
        distances = np.geomspace(max(h, 5.0), 300.0, 60)
        pga = [float(spectra.zhao_spectrum(6.5, r, h, source)["Sa_g"].iloc[0])
               for r in distances]
        ax.plot(distances, pga, color=COLORS[source], linewidth=2)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Distancia hipocentral (km)")
    ax.set_ylabel("PGA mediana (g), Mw 6.5")

    # Paso 4: curva de amenaza total y por fuente.
    ax = axes[3]
    pga_g = hazard["PGA_roca_cm_s2"] / 981.0
    for source in SOURCES:
        ax.plot(pga_g, hazard[f"lambda_{source}_por_anio"], color=COLORS[source],
                linewidth=1.5, linestyle="--")
    ax.plot(pga_g, hazard["lambda_total_por_anio"], color="black", linewidth=2,
            label="Total")
    ax.axhline(1 / 475, color="gray", linewidth=0.8, linestyle=":")
    ax.plot(params["PGA_g"], 1 / 475, "o", color="red", markersize=6)
    ax.annotate(f"Tr 475 años\nPGA = {params['PGA_g']:.3f} g",
                (params["PGA_g"], 1 / 475), xytext=(8, 8),
                textcoords="offset points", fontsize=9)
    ax.set_yscale("log")
    ax.set_xlim(0, 0.6)
    ax.set_ylim(1e-5, 10)
    ax.set_xlabel("PGA (g)")
    ax.set_ylabel("λ (PGA > y) (1/año)")

    # Paso 5: espectro de aceleración con meseta.
    ax = axes[4]
    t = np.linspace(0, 4, 801)
    sa = spectra.plateau_spectrum(t, params["PGA_g"], params["Sa_meseta_g"],
                                  params["T0_s"], params["Tc_s"])
    ax.plot(design["T_s"], design["Envolvente_respuestas_g"], color="gray",
            linestyle="--", linewidth=1.2)
    ax.plot(t, sa, color="red", linewidth=2.2)
    ax.set_xlim(0, 4)
    ax.set_ylim(bottom=0)
    ax.set_xlabel("Periodo T (s)")
    ax.set_ylabel("Sa (g)")

    handles = [plt.Line2D([], [], color=COLORS[s], linewidth=2.5, label=SHORT[s])
               for s in SOURCES]
    handles += [plt.Line2D([], [], color="black", linewidth=2, label="Amenaza total"),
                plt.Line2D([], [], color="red", linewidth=2.2,
                           label="Espectro de aceleración")]
    fig.legend(handles=handles, loc="lower center", ncol=5, frameon=False,
               fontsize=10, bbox_to_anchor=(0.5, -0.04))
    fig.tight_layout(rect=(0, 0.06, 1, 1), w_pad=1.5)
    output = FIGURE_DIRECTORY / "metodologia_cinco_pasos.png"
    fig.savefig(output, bbox_inches="tight")
    plt.close(fig)
    print(f"Figura: {output}")


if __name__ == "__main__":
    main()
