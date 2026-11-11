"""Espectros de respuesta en roca para Neiva (Paso 5).

Para cada sismofuente se toma el sismo de control de la desagregación de
475 años (Mw y distancia hipocentral medias) y se calcula el espectro de
respuesta mediano en roca con:

* Sadigh et al. (1997), roca y Mw <= 6.5, para la fuente cortical
  (ecuación desarrollada para sismos corticales someros).
* Zhao et al. (2006), Tablas 4 y 5, para la subducción interplaca e
  intraplaca (términos propios de cada tipo de fuente).

El espectro de diseño en roca es un espectro suavizado con meseta (rampa,
meseta y rama descendente 1/T, como en la NSR-10) que envuelve los espectros
de respuesta de cada fuente escalados al PGA de diseño de 475 años.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

DELIVERY_ROOT = Path(__file__).resolve().parents[1]
DATA_DIRECTORY = DELIVERY_ROOT / "datos" / "datos estadistica" / "integracion paso4"
DISAGGREGATION_PATH = DATA_DIRECTORY / "desagregacion_475_paso4.csv"
DESIGN_PATH = DATA_DIRECTORY / "pga_diseno_neiva_roca_periodos_retorno_paso4.csv"
OUTPUT_DIRECTORY = DELIVERY_ROOT / "datos" / "datos estadistica" / "espectros paso5"
FIGURE_DIRECTORY = DELIVERY_ROOT / "figuras"

GRAVITY_CM_S2 = 981.0
DEPTH_REFERENCE_KM = 15.0
DEPTH_CAP_KM = 125.0

# Zhao et al. (2006), Tabla 4 (a, b, c, d, e, S_I, S_S, S_SL) y Tabla 5
# (C1 roca, sigma_T). Periodo 0 = PGA. Unidades: cm/s2.
ZHAO_COLUMNS = ("T", "a", "b", "c", "d", "e", "SI", "SS", "SSL", "C1", "sigmaT")
ZHAO_TABLE = (
    (0.00, 1.101, -0.00564, 0.0055, 1.080, 0.01412, 0.000, 2.607, -0.528, 1.111, 0.723),
    (0.05, 1.076, -0.00671, 0.0075, 1.060, 0.01463, 0.000, 2.764, -0.551, 1.684, 0.779),
    (0.10, 1.118, -0.00787, 0.0090, 1.083, 0.01423, 0.000, 2.156, -0.420, 2.061, 0.849),
    (0.15, 1.134, -0.00722, 0.0100, 1.053, 0.01509, 0.000, 2.161, -0.431, 1.916, 0.839),
    (0.20, 1.147, -0.00659, 0.0120, 1.014, 0.01462, 0.000, 1.901, -0.372, 1.669, 0.811),
    (0.25, 1.149, -0.00590, 0.0140, 0.966, 0.01459, 0.000, 1.814, -0.360, 1.468, 0.786),
    (0.30, 1.163, -0.00520, 0.0150, 0.934, 0.01458, 0.000, 2.181, -0.450, 1.172, 0.770),
    (0.40, 1.200, -0.00422, 0.0100, 0.959, 0.01257, -0.041, 2.432, -0.506, 0.655, 0.766),
    (0.50, 1.250, -0.00338, 0.0060, 1.008, 0.01114, -0.053, 2.629, -0.554, 0.071, 0.760),
    (0.60, 1.293, -0.00282, 0.0030, 1.088, 0.01019, -0.103, 2.702, -0.575, -0.429, 0.766),
    (0.70, 1.336, -0.00258, 0.0025, 1.084, 0.00979, -0.146, 2.654, -0.572, -0.866, 0.769),
    (0.80, 1.386, -0.00242, 0.0022, 1.088, 0.00944, -0.164, 2.480, -0.540, -1.325, 0.770),
    (0.90, 1.433, -0.00232, 0.0020, 1.109, 0.00972, -0.206, 2.332, -0.522, -1.732, 0.771),
    (1.00, 1.479, -0.00220, 0.0020, 1.115, 0.01005, -0.239, 2.233, -0.509, -2.152, 0.775),
    (1.25, 1.551, -0.00207, 0.0020, 1.083, 0.01003, -0.256, 2.029, -0.469, -2.923, 0.773),
    (1.50, 1.621, -0.00224, 0.0020, 1.091, 0.00928, -0.306, 1.589, -0.379, -3.548, 0.779),
    (2.00, 1.694, -0.00201, 0.0025, 1.055, 0.00833, -0.321, 0.966, -0.248, -4.410, 0.787),
    (2.50, 1.748, -0.00187, 0.0028, 1.052, 0.00776, -0.337, 0.789, -0.221, -5.049, 0.786),
    (3.00, 1.759, -0.00147, 0.0032, 1.025, 0.00644, -0.331, 1.037, -0.263, -5.431, 0.776),
    (4.00, 1.826, -0.00195, 0.0040, 1.044, 0.00590, -0.390, 0.561, -0.169, -6.181, 0.751),
    (5.00, 1.825, -0.00237, 0.0050, 1.065, 0.00510, -0.498, 0.225, -0.120, -6.347, 0.745),
)

# Sadigh et al. (1997), roca, Mw <= 6.5 (tabla ASCE de clase). Unidades: g.
SADIGH_COLUMNS = ("T", "C1", "C2", "C3", "C4", "C5", "C6", "C7")
SADIGH_TABLE = (
    (0.00, -0.624, 1.0, 0.000, -2.100, 1.296, 0.250, 0.000),
    (0.07, 0.110, 1.0, 0.006, -2.128, 1.296, 0.250, -0.082),
    (0.10, 0.275, 1.0, 0.006, -2.148, 1.296, 0.250, -0.041),
    (0.20, 0.153, 1.0, -0.004, -2.080, 1.296, 0.250, 0.000),
    (0.30, -0.057, 1.0, -0.017, -2.028, 1.296, 0.250, 0.000),
    (0.40, -0.298, 1.0, -0.028, -1.990, 1.296, 0.250, 0.000),
    (0.50, -0.588, 1.0, -0.040, -1.945, 1.296, 0.250, 0.000),
    (0.75, -1.208, 1.0, -0.050, -1.865, 1.296, 0.250, 0.000),
    (1.00, -1.705, 1.0, -0.055, -1.800, 1.296, 0.250, 0.000),
    (1.50, -2.407, 1.0, -0.065, -1.725, 1.296, 0.250, 0.000),
    (2.00, -2.945, 1.0, -0.070, -1.670, 1.296, 0.250, 0.000),
    (3.00, -3.700, 1.0, -0.080, -1.610, 1.296, 0.250, 0.000),
    (4.00, -4.230, 1.0, -0.100, -1.570, 1.296, 0.250, 0.000),
)
SADIGH_MAX_MAGNITUDE = 6.5

SOURCE_TYPES = {
    "Cortical": "crustal",
    "Subducción interplaca": "interface",
    "Subducción intraplaca": "slab",
}
SOURCE_COLORS = {
    "Cortical": "#1f77b4",
    "Subducción interplaca": "#ff7f0e",
    "Subducción intraplaca": "#2ca02c",
}


def zhao_spectrum(magnitude: float, distance_km: float, depth_km: float,
                  source: str) -> pd.DataFrame:
    """Espectro mediano de Zhao et al. (2006) en roca (Clase I), en g."""
    table = pd.DataFrame(ZHAO_TABLE, columns=ZHAO_COLUMNS)
    kind = SOURCE_TYPES[source]
    depth = float(np.clip(depth_km, 0.0, DEPTH_CAP_KM))
    depth_term = table["e"] * (depth - DEPTH_REFERENCE_KM) if depth >= DEPTH_REFERENCE_KM else 0.0
    interface = table["SI"] if kind == "interface" else 0.0
    slab = table["SS"] + table["SSL"] * np.log(distance_km) if kind == "slab" else 0.0
    r = distance_km + table["c"] * np.exp(table["d"] * magnitude)
    ln_y = (
        table["a"] * magnitude
        + table["b"] * distance_km
        - np.log(r)
        + depth_term
        + interface
        + slab
        + table["C1"]
    )
    return pd.DataFrame({"T": table["T"], "Sa_g": np.exp(ln_y) / GRAVITY_CM_S2})


def sadigh_spectrum(magnitude: float, distance_km: float) -> pd.DataFrame:
    """Espectro mediano de Sadigh et al. (1997) en roca, Mw <= 6.5, en g.

    ln y = C1 + C2 M + C3 (8.5 - M)^2.5 + C4 ln(r + exp(C5 + C6 M))
           + C7 ln(r + 2)
    """
    if magnitude > SADIGH_MAX_MAGNITUDE:
        raise ValueError("La tabla de Sadigh usada solo aplica para Mw <= 6.5.")
    t = pd.DataFrame(SADIGH_TABLE, columns=SADIGH_COLUMNS)
    ln_y = (
        t["C1"]
        + t["C2"] * magnitude
        + t["C3"] * (8.5 - magnitude) ** 2.5
        + t["C4"] * np.log(distance_km + np.exp(t["C5"] + t["C6"] * magnitude))
        + t["C7"] * np.log(distance_km + 2.0)
    )
    return pd.DataFrame({"T": t["T"], "Sa_g": np.exp(ln_y)})


# Espectro NSR-10 de la Entrega 1: perfil B, Aa = Av = 0.25, Fa = Fv = 1, I = 1.
NSR10_AA = 0.25
NSR10_AV = 0.25
NSR10_PLATEAU_G = 2.5 * NSR10_AA


def nsr10_spectrum(periods: np.ndarray) -> np.ndarray:
    """Espectro elástico NSR-10 en roca usado como objetivo en la Entrega 1."""
    t = np.asarray(periods, dtype=float)
    t0, tc, tl = 0.1 * NSR10_AV / NSR10_AA, 0.48 * NSR10_AV / NSR10_AA, 2.4
    safe = np.maximum(t, 1e-9)
    return np.where(
        t < t0, NSR10_AA * (1 + 1.5 * t / t0),
        np.where(t <= tc, NSR10_PLATEAU_G,
                 np.where(t <= tl, 1.2 * NSR10_AV / safe, 1.2 * NSR10_AV * tl / safe**2)))


def plateau_spectrum(periods: np.ndarray, pga: float, sa_max: float,
                     t0: float, tc: float) -> np.ndarray:
    """Espectro de diseño con meseta: rampa lineal PGA -> Sa_max hasta T0,
    meseta Sa_max entre T0 y Tc, y rama Sa_max * Tc / T para T > Tc."""
    t = np.asarray(periods, dtype=float)
    ramp = pga + (sa_max - pga) * t / t0
    tail = sa_max * tc / np.maximum(t, 1e-9)
    return np.where(t < t0, ramp, np.where(t <= tc, sa_max, tail))


def fit_plateau(grid: np.ndarray, envelope: np.ndarray, pga: float) -> tuple[float, float, float]:
    """Ajusta (Sa_max, T0, Tc) para que el espectro con meseta quede por
    encima de la envolvente de respuestas en todos los periodos."""
    sa_max = float(np.nanmax(envelope))
    fine = np.linspace(0.0, grid.max(), 5001)
    env = np.interp(fine, grid, envelope)
    # T0: la rampa debe quedar por encima de la envolvente.
    rising = (fine > 0) & (env > pga)
    t0 = float(np.min((sa_max - pga) * fine[rising] / (env[rising] - pga)))
    t0 = np.floor(t0 * 100) / 100
    # Tc: Sa_max * Tc / T >= Sa(T) para todo T  ->  Tc = max(Sa * T) / Sa_max.
    tc = float(np.max(env * fine) / sa_max)
    tc = np.ceil(tc * 100) / 100
    return sa_max, t0, tc


def run_analysis() -> dict[str, pd.DataFrame]:
    scenarios = pd.read_csv(DISAGGREGATION_PATH, encoding="utf-8-sig")
    scenarios = scenarios[scenarios["Sismofuente"] != "Total"].set_index("Sismofuente")
    design = pd.read_csv(DESIGN_PATH, encoding="utf-8-sig")
    pga_design_g = float(
        design.loc[design["Periodo_retorno_anios"] == 475, "PGA_roca_g"].iloc[0]
    )

    spectra = {}
    rows = []
    for source in SOURCE_TYPES:
        s = scenarios.loc[source]
        m, r, h = float(s["Mw_media"]), float(s["R_media_km"]), float(s["profundidad_media_km"])
        if SOURCE_TYPES[source] == "crustal":
            spectra[(source, "Sadigh")] = sadigh_spectrum(m, r)
        else:
            spectra[(source, "Zhao")] = zhao_spectrum(m, r, h, source)
        rows.append({"Sismofuente": source, "Mw": m, "R_hipocentral_km": r, "h_km": h})
    scenario_table = pd.DataFrame(rows)

    # Espectro de diseño: formas espectrales escaladas al PGA de 475 años.
    grid = np.array([row[0] for row in ZHAO_TABLE])
    scaled = {}
    for key, spectrum in spectra.items():
        shape = spectrum["Sa_g"] / spectrum["Sa_g"].iloc[0]
        values = np.interp(grid, spectrum["T"], shape * pga_design_g)
        # Sin extrapolar fuera del rango de periodos de cada tabla.
        values[grid > spectrum["T"].max()] = np.nan
        scaled[key] = values
    # Envolvente: Sadigh para la fuente cortical y Zhao para la subducción.
    envelope = np.nanmax(np.vstack(list(scaled.values())), axis=0)

    # Espectro de diseño con meseta que envuelve las respuestas.
    sa_max, t0, tc = fit_plateau(grid, envelope, pga_design_g)
    design_values = plateau_spectrum(grid, pga_design_g, sa_max, t0, tc)
    # Aa equivalente: en la NSR-10 la meseta en roca vale 2.5 Aa.
    aa_equivalent = sa_max / 2.5

    OUTPUT_DIRECTORY.mkdir(parents=True, exist_ok=True)
    long_rows = []
    for (source, model), spectrum in spectra.items():
        for t, sa in zip(spectrum["T"], spectrum["Sa_g"]):
            long_rows.append({"Sismofuente": source, "Modelo": model, "T_s": t, "Sa_mediano_g": sa})
    response = pd.DataFrame(long_rows)
    response.to_csv(OUTPUT_DIRECTORY / "espectros_respuesta_fuentes_paso5.csv",
                    index=False, encoding="utf-8-sig")
    design_spectrum = pd.DataFrame({"T_s": grid, "Sa_diseno_meseta_g": design_values,
                                    "Envolvente_respuestas_g": envelope})
    for (source, model), values in scaled.items():
        design_spectrum[f"{source} ({model})"] = values
    design_spectrum.to_csv(OUTPUT_DIRECTORY / "espectro_diseno_roca_paso5.csv",
                           index=False, encoding="utf-8-sig")
    scenario_table.to_csv(OUTPUT_DIRECTORY / "escenarios_control_paso5.csv",
                          index=False, encoding="utf-8-sig")
    pd.DataFrame([{"PGA_g": pga_design_g, "Sa_meseta_g": sa_max, "T0_s": t0, "Tc_s": tc,
                   "Aa_equivalente": aa_equivalent}]).to_csv(
        OUTPUT_DIRECTORY / "parametros_espectro_diseno_paso5.csv",
        index=False, encoding="utf-8-sig")

    fine = np.linspace(0.0, 5.0, 2001)
    design_fine = plateau_spectrum(fine, pga_design_g, sa_max, t0, tc)
    nsr_fine = nsr10_spectrum(fine)
    plt.rcParams.update({"font.size": 12})

    # Figura 1 (punto 5): espectros de respuesta medianos de cada fuente.
    fig, ax = plt.subplots(figsize=(10, 6), dpi=300)
    for (source, model), spectrum in spectra.items():
        s = scenarios.loc[source]
        ax.plot(spectrum["T"], spectrum["Sa_g"], color=SOURCE_COLORS[source],
                linestyle="--" if model == "Sadigh" else "-",
                marker="o", markersize=5, linewidth=2.2,
                label=f"{source} – {model} (Mw {s['Mw_media']:.1f}, R {s['R_media_km']:.0f} km)")
    ax.set_xlim(0, 4)
    ax.set_ylim(bottom=0)
    ax.set_xlabel("Periodo T (s)")
    ax.set_ylabel("Sa mediana (g)")
    ax.set_title("Espectros de respuesta en roca por sismofuente — Neiva\n"
                 "Sismo de control de la desagregación de 475 años, 5 % de amortiguamiento",
                 fontsize=13)
    ax.grid(True, alpha=0.3)
    ax.legend(frameon=False, fontsize=11)
    fig.tight_layout()
    fig.savefig(FIGURE_DIRECTORY / "espectros_respuesta_fuentes_neiva.png", bbox_inches="tight")
    plt.close(fig)

    # Figura 2 (punto 6): espectro de aceleración (diseño) en roca.
    # Arriba: acercamiento 0-1 s. Abajo: rango completo 0-4 s.
    fig, axes = plt.subplots(2, 1, figsize=(10, 10), dpi=300,
                             gridspec_kw={"height_ratios": [1.15, 1]})
    y_top = max(sa_max, NSR10_PLATEAU_G) * 1.15
    for ax, x_max in zip(axes, (1.0, 4.0)):
        ax.axvspan(t0, tc, color="red", alpha=0.07, zorder=0)
        ax.plot(grid, envelope, color="gray", linestyle="--", marker="o", markersize=4,
                linewidth=1.8, label="Envolvente de los espectros de respuesta escalados")
        ax.plot(fine, design_fine, color="red", linewidth=3.2,
                label="Espectro de aceleración en roca (este estudio)")
        ax.plot(fine, nsr_fine, color="black", linestyle="-.", linewidth=2.0,
                label="NSR-10, perfil B (Entrega 1)")
        ax.set_xlim(0, x_max)
        ax.set_ylim(0, y_top)
        ax.set_xlabel("Periodo T (s)")
        ax.set_ylabel("Sa (g)")
        ax.grid(True, alpha=0.3)

    top, bottom = axes
    top.set_title("a) Acercamiento 0–1 s: rampa y meseta", loc="left", fontsize=12)
    top.annotate(f"PGA = {pga_design_g:.3f} g", xy=(0, pga_design_g),
                 xytext=(0.03, 0.06), fontsize=11,
                 arrowprops=dict(arrowstyle="->", color="black", lw=1))
    top.annotate(f"Meseta = {sa_max:.3f} g\nT0 = {t0:.2f} s, Tc = {tc:.2f} s",
                 xy=(tc, sa_max), xytext=(0.46, 0.315), fontsize=11, color="red",
                 arrowprops=dict(arrowstyle="->", color="red", lw=1))
    top.annotate(f"Meseta NSR-10 = {NSR10_PLATEAU_G:.3f} g",
                 xy=(0.40, NSR10_PLATEAU_G), xytext=(0.55, NSR10_PLATEAU_G + 0.04),
                 fontsize=11, arrowprops=dict(arrowstyle="->", color="black", lw=1))
    bottom.set_title("b) Rango completo 0–4 s", loc="left", fontsize=12)

    handles, labels = top.get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=1, frameon=False,
               fontsize=11, bbox_to_anchor=(0.5, -0.01))
    fig.suptitle("Espectro de aceleración en roca — Neiva\n"
                 "Tr = 475 años, 5 % de amortiguamiento", fontsize=14)
    fig.tight_layout(rect=(0, 0.09, 1, 1))
    fig.savefig(FIGURE_DIRECTORY / "espectro_diseno_roca_neiva.png", bbox_inches="tight")
    plt.close(fig)
    plt.rcParams.update(plt.rcParamsDefault)

    # Verificación: la meseta debe quedar por encima de todas las respuestas.
    check = np.interp(fine, grid, envelope)
    margin = plateau_spectrum(fine, pga_design_g, sa_max, t0, tc) - check
    assert margin.min() > -1e-9, f"El espectro de diseño no envuelve: {margin.min()}"

    print("Escenarios de control:")
    print(scenario_table.to_string(index=False))
    print("\nEspectros medianos (g):")
    print(response.pivot_table(index="T_s", columns=["Sismofuente", "Modelo"],
                               values="Sa_mediano_g").round(4).to_string())
    print("\nEspectro de diseño (g):")
    print(design_spectrum.round(3).to_string(index=False))
    print(f"\nMeseta: Sa_max = {sa_max:.3f} g, T0 = {t0:.2f} s, Tc = {tc:.2f} s")
    print(f"Aa equivalente = Sa_max / 2.5 = {aa_equivalent:.3f}")
    return {"respuesta": response, "diseno": design_spectrum}


if __name__ == "__main__":
    run_analysis()
