"""Curvas normalizadas G/Gmax y amortiguamiento por estrato (punto 2).

E1 (arena): ajuste a los ensayos. La columna resonante de 10 m da Gmax y la
parte de deformaciones pequeñas; el triaxial cíclico de 5 m, la de
deformaciones grandes (ambos con sigma3 = 200 kPa).

E2 (arcilla CL, CH en la base): forma de Darendeli (2001), con a = 0.919,
ajustada a la columna resonante de 20 m (solo llega a 0.05 %; la forma
extrapola a deformaciones grandes).

La figura muestra, como en la tarea de referencia, G/Gmax y D en la misma
gráfica (dos ejes) con los rangos de deformación de cada ensayo.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

DELIVERY_ROOT = Path(__file__).resolve().parents[1]
LAB = DELIVERY_ROOT / "datos" / "ensayos_laboratorio"
RESULTS = DELIVERY_ROOT / "datos" / "resultados"
FIGURES = DELIVERY_ROOT / "figuras"

GAMMA_W = 9.81          # kN/m3; nivel freático en superficie
K0 = 0.5                # coeficiente de empuje en reposo supuesto
P_ATM = 101.3           # kPa
N_CYCLES = 10
FREQUENCY_HZ = 1.0
STRAINS = np.logspace(-5, 1, 61)   # deformación cortante en %

# Estratos: (nombre, techo, base, peso unitario, IP representativo)
STRATA = {
    "E1": (0.0, 16.5, 18.6, 5.5),    # IP de M3 y M6 (6 y 5)
    "E2": (16.5, 30.0, 18.7, 17.8),  # IP de M16, M22, M26 y M30 (11, 16, 17, 27)
}


def mean_effective_stress(depth: float) -> float:
    """Esfuerzo medio efectivo (kPa) a una profundidad, con NF en superficie."""
    sigma_v = 0.0
    for top, bottom, gamma, _ in STRATA.values():
        h = max(0.0, min(depth, bottom) - top)
        sigma_v += h * (gamma - GAMMA_W)
    return sigma_v * (1 + 2 * K0) / 3


def hyperbolic(gamma, gamma_r, a):
    return 1.0 / (1.0 + (gamma / gamma_r) ** a)


def masing_damping(gamma, gamma_r, a):
    """Amortiguamiento de Masing (%) según Darendeli (2001)."""
    d1 = (100 / np.pi) * (4 * (gamma - gamma_r * np.log((gamma + gamma_r) / gamma_r))
                          / (gamma**2 / (gamma + gamma_r)) - 2)
    c1 = -1.1143 * a**2 + 1.8618 * a + 0.2523
    c2 = 0.0805 * a**2 - 0.0710 * a - 0.0095
    c3 = -0.0005 * a**2 + 0.0002 * a + 0.0003
    return c1 * d1 + c2 * d1**2 + c3 * d1**3


def darendeli(gamma, pi, sigma_m_kpa, ocr=1.0):
    """Curvas de Darendeli (2001): G/Gmax y D (%) contra gamma (%)."""
    s = sigma_m_kpa / P_ATM
    gamma_r = (0.0352 + 0.0010 * pi * ocr**0.3246) * s**0.3483
    a = 0.919
    d_min = (0.8005 + 0.0129 * pi * ocr**-0.1069) * s**-0.2889 * (1 + 0.2919 * np.log(FREQUENCY_HZ))
    b = 0.6329 - 0.0057 * np.log(N_CYCLES)
    g_ratio = hyperbolic(gamma, gamma_r, a)
    damping = b * g_ratio**0.1 * masing_damping(gamma, gamma_r, a) + d_min
    return g_ratio, damping, gamma_r, d_min


def run_analysis() -> dict[str, pd.DataFrame]:
    rc = pd.read_csv(LAB / "columna_resonante.csv", encoding="utf-8-sig")
    tx = pd.read_csv(LAB / "triaxial_ciclico_5m.csv", encoding="utf-8-sig")

    # Gmax de campo de cada estrato (down hole): rho del estrato y Vs equivalente
    # por tiempo de viaje, Vs = sum(h) / sum(h / Vs_capa).
    tramos = pd.read_csv(RESULTS / "perfil_tramos_vs_gmax.csv", encoding="utf-8-sig")
    gmax_dh = {}
    for name, sub in tramos.groupby("estrato"):
        h = sub["hasta_m"] - sub["desde_m"]
        vs_eq = h.sum() / (h / sub["Vs_m_s"]).sum()
        rho = sub["peso_unitario_kN_m3"].iloc[0] / 9.81 * 1000
        gmax_dh[name] = rho * vs_eq**2 / 1e6   # MPa

    # --- E1: datos de laboratorio normalizados con el Gmax del down hole.
    rc10 = rc[rc["profundidad_m"] == 10.0]
    g_max_e1 = gmax_dh["E1"] * 1000.0  # kPa
    lab_e1 = pd.concat([
        pd.DataFrame({"ensayo": "Columna resonante 10 m", "gamma_pct": rc10["gamma_pct"],
                      "G_Gmax": rc10["G_MPa"] * 1000 / g_max_e1, "D_pct": rc10["D_pct"]}),
        pd.DataFrame({"ensayo": "Triaxial cíclico 5 m", "gamma_pct": tx["gamma_pct"],
                      "G_Gmax": tx["G_kPa"] / g_max_e1, "D_pct": tx["D_pct"]}),
    ], ignore_index=True)
    # Forma de la curva: G dividido por el módulo inicial del ensayo (columna
    # resonante a 0.002 %); el Gmax de campo es el del down hole.
    g0_lab_e1 = float(rc10["G_MPa"].iloc[0]) * 1000
    lab_e1["G_G0"] = lab_e1["G_Gmax"] * g_max_e1 / g0_lab_e1
    (gr_e1, a_e1), _ = curve_fit(hyperbolic, lab_e1["gamma_pct"], lab_e1["G_G0"],
                                 p0=(0.1, 1.0), bounds=([1e-4, 0.3], [10, 2]))
    g_lab = hyperbolic(lab_e1["gamma_pct"].to_numpy(), gr_e1, a_e1)

    def damping_model(g_ratio, d_min, d_max, n):
        return d_min + (d_max - d_min) * (1 - g_ratio) ** n

    (dmin_e1, dmax_e1, n_e1), _ = curve_fit(
        damping_model, g_lab, lab_e1["D_pct"], p0=(0.5, 18, 1),
        bounds=([0, 5, 0.3], [3, 40, 3]))
    r2_g = 1 - np.sum((lab_e1["G_G0"] - g_lab) ** 2) / np.sum(
        (lab_e1["G_G0"] - lab_e1["G_G0"].mean()) ** 2)

    curves = {"gamma_pct": STRAINS}
    g_e1 = hyperbolic(STRAINS, gr_e1, a_e1)
    curves["E1_G_Gmax"] = g_e1
    curves["E1_D_pct"] = damping_model(g_e1, dmin_e1, dmax_e1, n_e1)

    params = [{
        "estrato": "E1", "metodo": "Ajuste a ensayos (CR 10 m + TC 5 m)",
        "IP": STRATA["E1"][3], "sigma_m_kPa": 200.0, "gamma_r_pct": gr_e1, "a": a_e1,
        "D_min_pct": dmin_e1, "R2_G": r2_g,
    }]

    # --- E2: forma de Darendeli (a = 0.919) ajustada a la columna resonante de 20 m.
    # La muestra da un módulo inicial (264 MPa) menor que el Gmax de campo del down
    # hole (286 MPa), como es usual por la alteración de la muestra. La forma de la
    # curva se toma del ensayo, normalizado con su propio módulo inicial G0, y el
    # Gmax es el del down hole. En la figura los puntos se muestran con G/Gmax_DH.
    rc20 = rc[rc["profundidad_m"] == 20.0].copy()
    rc20["G_Gmax"] = rc20["G_MPa"] / gmax_dh["E2"]
    rc20["G_G0"] = rc20["G_MPa"] / rc20["G_MPa"].iloc[0]
    a_dar = 0.919
    (gr_e2,), _ = curve_fit(lambda g, gr: hyperbolic(g, gr, a_dar), rc20["gamma_pct"],
                            rc20["G_G0"], p0=(0.1,), bounds=([1e-4], [10]))
    b_dar = 0.6329 - 0.0057 * np.log(N_CYCLES)

    def clay_damping(gamma, d_min):
        g = hyperbolic(gamma, gr_e2, a_dar)
        return b_dar * g**0.1 * masing_damping(gamma, gr_e2, a_dar) + d_min

    (dmin_e2,), _ = curve_fit(clay_damping, rc20["gamma_pct"], rc20["D_pct"],
                              p0=(0.1,), bounds=([0.0], [5.0]))
    g_e2 = hyperbolic(STRAINS, gr_e2, a_dar)
    curves["E2_G_Gmax"] = g_e2
    curves["E2_D_pct"] = clay_damping(STRAINS, dmin_e2)
    params.append({"estrato": "E2", "metodo": "Forma de Darendeli ajustada a CR 20 m",
                   "IP": STRATA["E2"][3], "sigma_m_kPa": 400.0, "gamma_r_pct": gr_e2,
                   "a": a_dar, "D_min_pct": dmin_e2, "R2_G": np.nan})

    # Comparación: Darendeli puro con las condiciones de cada ensayo.
    g_check, d_check, gr_check, _ = darendeli(rc20["gamma_pct"].to_numpy(), STRATA["E2"][3], 400.0)
    rc20["G_Gmax_Darendeli_400kPa"] = g_check
    rc20["D_Darendeli_400kPa"] = d_check
    _, _, gr_e1_dar, _ = darendeli(np.array([0.01]), STRATA["E1"][3], 200.0)
    print(f"gamma_r Darendeli: E1 (200 kPa) = {gr_e1_dar:.3f} %, E2 (400 kPa) = {gr_check:.3f} %")

    curves_df = pd.DataFrame(curves)
    params_df = pd.DataFrame(params)
    RESULTS.mkdir(parents=True, exist_ok=True)
    curves_df.to_csv(RESULTS / "curvas_degradacion_por_estrato.csv", index=False, encoding="utf-8-sig")
    params_df.to_csv(RESULTS / "parametros_curvas_degradacion.csv", index=False, encoding="utf-8-sig")
    lab_e1.to_csv(RESULTS / "puntos_laboratorio_E1.csv", index=False, encoding="utf-8-sig")
    rc20.to_csv(RESULTS / "verificacion_columna_resonante_20m.csv", index=False, encoding="utf-8-sig")

    # --- Figura: G/Gmax (eje izquierdo) y D (eje derecho) por estrato.
    plt.rcParams.update({"font.size": 11})
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5), dpi=300)
    blue, red = "#1f4e79", "#c0392b"
    ranges = (("Down hole", 1e-5, 1e-4), ("Columna resonante", 1e-4, 1e-2),
              ("Triaxial cíclico", 1e-2, 10.0))
    tests = {
        "E1": [(lab_e1[lab_e1["ensayo"] == t], m, t) for t, m in
               (("Columna resonante 10 m", "o"), ("Triaxial cíclico 5 m", "s"))],
        "E2": [(rc20, "^", "Columna resonante 20 m")],
    }
    titles = {"E1": "E1: arena con grava y limo (SP-SM)",
              "E2": "E2: arcilla CL (CH en la base)"}
    for ax, name in zip(axes, ("E1", "E2")):
        axd = ax.twinx()
        for i, (txt, x0, x1) in enumerate(ranges):
            ax.axvspan(x0, x1, color=("#f2f2f2", "#e6e6e6")[i % 2], zorder=0)
            ax.text(np.sqrt(x0 * x1), 1.08, txt, ha="center", fontsize=8.5)
        # Tramo con datos (línea continua) y tramo extrapolado con el modelo (punteado).
        g_lim = 10.0 if name == "E1" else float(rc20["gamma_pct"].max())
        ok = STRAINS <= g_lim
        ax.plot(STRAINS[ok], curves[f"{name}_G_Gmax"][ok], color=blue, lw=2.2, label="G/Gmax")
        axd.plot(STRAINS[ok], curves[f"{name}_D_pct"][ok], color=red, lw=2.2, label="D")
        if (~ok).any():
            ext = STRAINS >= g_lim
            ax.plot(STRAINS[ext], curves[f"{name}_G_Gmax"][ext], color=blue, lw=2.2, ls="--")
            axd.plot(STRAINS[ext], curves[f"{name}_D_pct"][ext], color=red, lw=2.2, ls="--")
            ax.plot([], [], color="gray", lw=2, ls="--", label="Extrapolado (Darendeli)")
        ax.plot(3e-5, 1.0, "D", color=blue, ms=7,
                label=f"Down hole: Gmax = {gmax_dh[name]:.1f} MPa")
        for sub, mk, lab in tests[name]:
            ax.plot(sub["gamma_pct"], sub["G_G0"], mk, color=blue, mfc="white", mew=1.6, ms=7,
                    label=lab)
            axd.plot(sub["gamma_pct"], sub["D_pct"], mk, color=red, mfc="white", mew=1.6, ms=7)
        ax.set_xscale("log")
        ax.set_xlim(1e-5, 10)
        ax.set_ylim(0, 1.15)
        ax.set_yticks(np.arange(0, 1.01, 0.2))
        axd.set_ylim(0, 25)
        ax.set_xlabel("Deformación cortante γ (%)")
        ax.set_ylabel("G/Gmax", color=blue)
        axd.set_ylabel("Amortiguamiento D (%)", color=red)
        ax.tick_params(axis="y", colors=blue)
        axd.tick_params(axis="y", colors=red)
        ax.grid(True, which="major", alpha=0.3)
        ax.set_title(titles[name], fontsize=11, pad=22)
        h1, l1 = ax.get_legend_handles_labels()
        h2, l2 = axd.get_legend_handles_labels()
        ax.legend(h1 + h2, l1 + l2, frameon=False, fontsize=8.5, loc="center left")
    fig.tight_layout()
    fig.savefig(FIGURES / "curvas_degradacion_estratos.png", bbox_inches="tight")
    plt.close(fig)

    print(params_df.round(4).to_string(index=False))
    print("\nVerificación columna resonante 20 m:")
    print(rc20[["gamma_pct", "G_Gmax", "G_Gmax_Darendeli_400kPa", "D_pct",
                "D_Darendeli_400kPa"]].round(3).to_string(index=False))
    print("\nE1 laboratorio:")
    print(lab_e1.round(3).to_string(index=False))
    print(f"\nE1: Gmax = {g_max_e1/1000:.1f} MPa, D = {dmin_e1:.2f} + ({dmax_e1:.1f} - Dmin)(1-G/Gmax)^{n_e1:.2f}")
    return {"curvas": curves_df, "parametros": params_df}


if __name__ == "__main__":
    run_analysis()
