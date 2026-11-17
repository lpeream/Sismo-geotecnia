"""Perfil estratigráfico del sondeo P-1 (punto 1 del enunciado).

Usa el perfil de Vs obtenido al procesar el down hole
(datos/resultados/downhole_perfil_velocidades.csv, generado por
procesamiento_downhole.py) y la clasificación de las muestras. Calcula Vs30 y
Gmax = rho * Vs^2 por tramo y genera la figura del perfil, con la gráfica de
referencia del enunciado como comparación.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

DELIVERY_ROOT = Path(__file__).resolve().parents[1]
LAB = DELIVERY_ROOT / "datos" / "ensayos_laboratorio"
RESULTS = DELIVERY_ROOT / "datos" / "resultados"
FIGURES = DELIVERY_ROOT / "figuras"

# Estratos propuestos a partir de la clasificación y del down hole.
STRATA = (
    ("E1", 0.0, 16.5, "Arena con grava y limo (SP-SM)", "#e8d9a8"),
    ("E2", 16.5, 30.0, "Arcilla CL\n(CH en la base)", "#b9c9a3"),
)


def run_analysis() -> dict[str, pd.DataFrame]:
    vs = pd.read_csv(RESULTS / "downhole_perfil_velocidades.csv", encoding="utf-8-sig")
    ref = pd.read_csv(LAB / "vs_downhole_leido_de_grafica.csv", encoding="utf-8-sig")
    samples = pd.read_csv(LAB / "clasificacion_muestras.csv", encoding="utf-8-sig")

    travel_time = ((vs["hasta_m"] - vs["desde_m"]) / vs["Vs_m_s"]).sum()
    vs30 = 30.0 / travel_time

    # Peso unitario por estrato: promedio de las muestras de cada uno.
    # E1: triaxial CU 4 m (18.2), triaxial cíclico 5 m (18.6), CR 10 m (18.9).
    # E2: CR 20 m (19.3) y triaxial CU 28 m (18.1).
    unit_weight = {"E1": 18.6, "E2": 18.7}
    rows = []
    for _, r in vs.iterrows():
        # Divide el tramo de Vs en los límites de estrato que caen dentro de él.
        cuts = sorted({r["desde_m"], r["hasta_m"]}
                      | {s[1] for s in STRATA if r["desde_m"] < s[1] < r["hasta_m"]})
        for top, bottom in zip(cuts[:-1], cuts[1:]):
            mid = (top + bottom) / 2
            stratum = next(s for s in STRATA if s[1] <= mid < s[2])
            gamma = unit_weight[stratum[0]]
            rows.append({
                "desde_m": top, "hasta_m": bottom, "estrato": stratum[0],
                "Vs_m_s": r["Vs_m_s"], "peso_unitario_kN_m3": gamma,
                "Gmax_MPa": gamma / 9.81 * r["Vs_m_s"] ** 2 / 1000.0,
            })
    profile = pd.DataFrame(rows)
    RESULTS.mkdir(parents=True, exist_ok=True)
    profile.to_csv(RESULTS / "perfil_tramos_vs_gmax.csv", index=False, encoding="utf-8-sig")

    plt.rcParams.update({"font.size": 11})
    fig, (ax_vs, ax_soil) = plt.subplots(1, 2, figsize=(9, 7), dpi=300, sharey=True,
                                         gridspec_kw={"width_ratios": [1.3, 1]})
    for name, top, bottom, label, color in STRATA:
        for ax in (ax_vs, ax_soil):
            ax.axhspan(top, bottom, color=color, alpha=0.55, zorder=0)
        ax_soil.text(0.97, (top + bottom) / 2, f"{name}\n{label}", ha="right",
                     va="center", fontsize=9)
    depth = [v for _, r in vs.iterrows() for v in (r["desde_m"], r["hasta_m"])]
    speed = [v for _, r in vs.iterrows() for v in (r["Vs_m_s"], r["Vs_m_s"])]
    ax_vs.plot(speed, depth, color="black", linewidth=2, label="Este trabajo")
    ref_depth = [v for _, r in ref.iterrows() for v in (r["desde_m"], r["hasta_m"])]
    ref_speed = [v for _, r in ref.iterrows() for v in (r["Vs_m_s"], r["Vs_m_s"])]
    ax_vs.plot(ref_speed, ref_depth, color="gray", linewidth=1.6, linestyle="--",
               label="Referencia del enunciado")
    ax_vs.legend(frameon=False, loc="lower left", fontsize=9)
    ax_vs.set_xlabel("Vs (m/s)")
    ax_vs.set_ylabel("Profundidad (m)")
    ax_vs.set_xlim(0, 450)
    ax_vs.grid(True, alpha=0.3)
    ax_vs.set_title(f"Down hole (Vs30 ≈ {vs30:.0f} m/s)", fontsize=11)
    for _, s in samples.iterrows():
        ax_soil.plot(0.08, s["profundidad_m"], "o", color="black", markersize=4)
        ax_soil.annotate(f"{s['muestra']}: {s['USCS_informado']}", (0.08, s["profundidad_m"]),
                         xytext=(6, 0), textcoords="offset points", va="center", fontsize=8)
    ax_soil.set_xlim(0, 1)
    ax_soil.set_xticks([])
    ax_soil.set_title("Estratos y muestras", fontsize=11)
    ax_vs.set_ylim(30.6, 0)
    fig.tight_layout()
    fig.savefig(FIGURES / "perfil_estratigrafico_p1.png", bbox_inches="tight")
    plt.close(fig)

    print(f"Vs30 = {vs30:.1f} m/s")
    print(profile.round(1).to_string(index=False))
    return {"perfil": profile}


if __name__ == "__main__":
    run_analysis()
