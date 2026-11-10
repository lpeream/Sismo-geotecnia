from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.special import ndtr


ROOT = Path(__file__).resolve().parents[2]
ENTREGA_2 = ROOT / "Entrega_02"
CATALOG_PATH = (
    ENTREGA_2
    / "recursos"
    / "Qgis"
    / "datos"
    / "Sismos_260_clasificados_20261004_192107.csv"
)
ATTENUATION_PATH = (
    ENTREGA_2
    / "datos"
    / "datos estadistica"
    / "parametros_atenuacion_neiva_pga.json"
)
DATA_DIRECTORY = (
    ENTREGA_2
    / "datos"
    / "datos estadistica"
    / "integracion paso4"
)
FIGURE_DIRECTORY = ENTREGA_2 / "figuras"

RECURRENCE: dict[str, dict[str, float]] = {
    "Cortical": {"a": 4.197, "b": 0.882, "Mw_max": 7.2},
    "Subducción interplaca": {"a": 3.911, "b": 0.964, "Mw_max": 8.0},
    "Subducción intraplaca": {"a": 4.025, "b": 0.821, "Mw_max": 7.8},
}
SOURCE_TYPES = {
    "Cortical": "crustal",
    "Subducción interplaca": "interface",
    "Subducción intraplaca": "slab",
}
MINIMUM_MAGNITUDE = 3.0
MAGNITUDE_STEP = 0.025
SIGMA_LN = 0.65
PGA_DISPLAY_MAX_CM_S2 = 100.0
PGA_HAZARD_FIGURE_MAX_CM_S2 = 600.0
PGA_CALCULATION_MAX_CM_S2 = 1200.0
PGA_STEP_CM_S2 = 0.25
PGA_EXTENDED_STEP_CM_S2 = 0.5
RETURN_PERIODS_YEARS = (100, 475, 975, 2475)
EARTH_RADIUS_KM = 6371.0
MINIMUM_DISTANCE_KM = 0.1


def _read_inputs() -> tuple[pd.DataFrame, dict[str, Any], dict[str, Any]]:
    if not CATALOG_PATH.is_file():
        raise FileNotFoundError(f"No se encontró el catálogo: {CATALOG_PATH}")
    if not ATTENUATION_PATH.is_file():
        raise FileNotFoundError(
            f"No se encontró la configuración de atenuación: {ATTENUATION_PATH}"
        )

    catalog = pd.read_csv(CATALOG_PATH, encoding="latin1")
    metadata = json.loads(ATTENUATION_PATH.read_text(encoding="utf-8"))
    configuration = metadata["configuration"]
    parameters = metadata["model_parameters"]
    if configuration["site_class"] != "I":
        raise ValueError("La integración requiere atenuación en roca, Clase I.")
    if float(parameters["site_terms"]["I"]) != 0.0:
        raise ValueError("La Clase I debe usar el término de sitio C1 = 0.0.")
    if parameters["label"] != "PGA":
        raise ValueError("La configuración de atenuación debe corresponder a PGA.")

    required = {"latitud", "longitud", "Z", "Mw", "fuente"}
    missing = required.difference(catalog.columns)
    if missing:
        raise ValueError(
            "Faltan columnas requeridas en el catálogo: "
            + ", ".join(sorted(missing))
        )
    for column in ("latitud", "longitud", "Z", "Mw"):
        catalog[column] = pd.to_numeric(catalog[column], errors="coerce")
    if catalog[list(required)].isna().any().any():
        raise ValueError("El catálogo contiene datos no válidos para la integración.")
    unexpected_sources = set(catalog["fuente"]) - set(RECURRENCE)
    if unexpected_sources:
        raise ValueError(
            f"El catálogo contiene fuentes no configuradas: {sorted(unexpected_sources)}"
        )
    if (catalog["Mw"] < MINIMUM_MAGNITUDE).any():
        raise ValueError(
            "El catálogo de integración debe tener magnitudes Mw >= 3.0."
        )
    return catalog, configuration, parameters


def _epicentral_distance_km(
    latitude: np.ndarray,
    longitude: np.ndarray,
    site_latitude: float,
    site_longitude: float,
) -> np.ndarray:
    latitude_rad = np.radians(latitude)
    longitude_rad = np.radians(longitude)
    site_latitude_rad = np.radians(site_latitude)
    site_longitude_rad = np.radians(site_longitude)
    delta_latitude = latitude_rad - site_latitude_rad
    delta_longitude = longitude_rad - site_longitude_rad
    haversine = (
        np.sin(delta_latitude / 2.0) ** 2
        + np.cos(site_latitude_rad)
        * np.cos(latitude_rad)
        * np.sin(delta_longitude / 2.0) ** 2
    )
    central_angle = 2.0 * np.arcsin(np.sqrt(np.clip(haversine, 0.0, 1.0)))
    return EARTH_RADIUS_KM * central_angle


def _attenuation_median(
    magnitude: float,
    distances_km: np.ndarray,
    depth_km: float,
    source: str,
    parameters: dict[str, Any],
    depth_reference_km: float,
    depth_cap_km: float,
) -> np.ndarray:
    distance = np.maximum(distances_km, MINIMUM_DISTANCE_KM)
    depth = float(np.clip(depth_km, 0.0, depth_cap_km))
    tectonic_type = SOURCE_TYPES[source]
    interface_term = (
        float(parameters["interface"]) if tectonic_type == "interface" else 0.0
    )
    slab_term = float(parameters["slab"]) if tectonic_type == "slab" else 0.0
    slab_distance_term = (
        float(parameters["slab_distance"])
        if tectonic_type == "slab"
        else 0.0
    )
    depth_term = (
        float(parameters["e"]) * (depth - depth_reference_km)
        if depth >= depth_reference_km
        else 0.0
    )
    adjusted_distance = distance + float(parameters["c"]) * np.exp(
        float(parameters["d"]) * magnitude
    )
    ln_median = (
        float(parameters["a"]) * magnitude
        + float(parameters["b"]) * distance
        - np.log(adjusted_distance)
        + depth_term
        + interface_term
        + slab_term
        + slab_distance_term * np.log(distance)
        + float(parameters["site_terms"]["I"])
    )
    return np.exp(ln_median)


def _source_hazard(
    source: str,
    source_events: pd.DataFrame,
    pga: np.ndarray,
    parameters: dict[str, Any],
    configuration: dict[str, Any],
) -> tuple[np.ndarray, dict[str, float]]:
    recurrence = RECURRENCE[source]
    maximum_magnitude = recurrence["Mw_max"]
    if maximum_magnitude <= MINIMUM_MAGNITUDE:
        raise ValueError(f"Mw_max debe superar 3.0 para {source}.")
    if source_events.empty:
        raise ValueError(f"No hay registros espaciales para la fuente {source}.")

    edges = np.arange(
        MINIMUM_MAGNITUDE,
        maximum_magnitude,
        MAGNITUDE_STEP,
        dtype=float,
    )
    if edges.size == 0 or not np.isclose(edges[0], MINIMUM_MAGNITUDE):
        edges = np.array([MINIMUM_MAGNITUDE], dtype=float)
    if np.isclose(edges[-1], maximum_magnitude):
        edges[-1] = maximum_magnitude
    else:
        edges = np.append(edges, maximum_magnitude)
    lower = edges[:-1]
    upper = edges[1:]
    magnitude_midpoints = (lower + upper) / 2.0
    annual_rate_in_bins = 10.0 ** (
        recurrence["a"] - recurrence["b"] * lower
    ) - 10.0 ** (recurrence["a"] - recurrence["b"] * upper)
    minimum_rate = 10.0 ** (
        recurrence["a"] - recurrence["b"] * MINIMUM_MAGNITUDE
    )
    maximum_rate = 10.0 ** (
        recurrence["a"] - recurrence["b"] * maximum_magnitude
    )
    expected_truncated_rate = minimum_rate - maximum_rate
    if not np.isclose(
        annual_rate_in_bins.sum(),
        expected_truncated_rate,
        rtol=1e-10,
        atol=1e-12,
    ):
        raise ArithmeticError(
            f"Las tasas incrementales no conservan la tasa truncada de {source}."
        )

    distances = source_events["distancia_epicentral_km"].to_numpy(dtype=float)
    reference_depth = float(source_events["Z"].clip(lower=0).median())
    reference_depth = float(
        np.clip(reference_depth, 0.0, configuration["depth_cap_km"])
    )
    source_curve = np.zeros_like(pga)
    positive_pga = pga > 0.0
    log_pga = np.log(pga[positive_pga])
    for magnitude, annual_rate in zip(magnitude_midpoints, annual_rate_in_bins):
        median = _attenuation_median(
            float(magnitude),
            distances,
            reference_depth,
            source,
            parameters,
            float(configuration["depth_reference_km"]),
            float(configuration["depth_cap_km"]),
        )
        standardized = (
            log_pga[:, np.newaxis] - np.log(median)[np.newaxis, :]
        ) / SIGMA_LN
        conditional_exceedance = ndtr(-standardized).mean(axis=1)
        source_curve[positive_pga] += annual_rate * conditional_exceedance
    source_curve[0] = expected_truncated_rate

    geometry = {
        "N_catalogo": int(len(source_events)),
        "Mw_max": maximum_magnitude,
        "profundidad_mediana_km": reference_depth,
        "distancia_min_km": float(distances.min()),
        "distancia_p10_km": float(np.quantile(distances, 0.10)),
        "distancia_mediana_km": float(np.median(distances)),
        "distancia_p90_km": float(np.quantile(distances, 0.90)),
        "distancia_max_km": float(distances.max()),
        "tasa_truncada_Mw_ge_3_por_anio": expected_truncated_rate,
    }
    return source_curve, geometry


def _fit_exponential(pga: np.ndarray, rate: np.ndarray) -> dict[str, float]:
    if np.any(rate <= 0.0):
        raise ValueError("El ajuste exponencial requiere tasas estrictamente positivas.")
    slope, intercept = np.polyfit(pga, np.log(rate), 1)
    if slope >= 0.0:
        raise ValueError("La tasa debe decrecer con PGA para el ajuste exponencial.")
    log_rate = np.log(rate)
    fitted_log_rate = intercept + slope * pga
    total_squares = float(np.sum((log_rate - log_rate.mean()) ** 2))
    residual_squares = float(np.sum((log_rate - fitted_log_rate) ** 2))
    r_squared = 1.0 - residual_squares / total_squares
    return {
        "A_eventos_por_anio": float(np.exp(intercept)),
        "B_por_cm_s2": float(-slope),
        "R2_ln_lambda": r_squared,
    }


def _design_accelerations(
    pga: np.ndarray, total_rate: np.ndarray
) -> pd.DataFrame:
    rows = []
    for return_period in RETURN_PERIODS_YEARS:
        target_rate = 1.0 / return_period
        if target_rate > total_rate[0] or target_rate < total_rate[-1]:
            raise ValueError(
                f"La tasa objetivo para Tr={return_period} años no queda "
                "dentro del rango calculado para PGA."
            )
        design_pga = float(
            np.interp(
                np.log(target_rate),
                np.log(total_rate[::-1]),
                pga[::-1],
            )
        )
        rows.append(
            {
                "Periodo_retorno_anios": return_period,
                "Tasa_objetivo_anual_1_Tr": target_rate,
                "PGA_roca_cm_s2": design_pga,
                "PGA_roca_g": design_pga / 981.0,
                "P_excedencia_50_anios": 1.0
                - np.exp(-50.0 * target_rate),
            }
        )
    return pd.DataFrame(rows)


def _plot_hazard(
    pga: np.ndarray,
    source_curves: dict[str, np.ndarray],
    total_curve: np.ndarray,
    design: pd.DataFrame,
) -> Path:
    colors = {
        "Cortical": "#1f77b4",
        "Subducción interplaca": "#ff7f0e",
        "Subducción intraplaca": "#2ca02c",
    }
    figure, axis = plt.subplots(figsize=(10, 6.2), dpi=300)
    visible = pga <= PGA_HAZARD_FIGURE_MAX_CM_S2
    for source, curve in source_curves.items():
        axis.plot(
            pga[visible],
            curve[visible],
            linestyle=":",
            linewidth=1.8,
            color=colors[source],
            label=source,
        )
    axis.plot(
        pga[visible],
        total_curve[visible],
        color="black",
        linewidth=2.5,
        label="Amenaza total en roca",
    )
    design_475 = design.loc[design["Periodo_retorno_anios"] == 475]
    if len(design_475) != 1:
        raise ValueError("Se requiere un único resultado de diseño para Tr=475 años.")
    design_pga = float(design_475.iloc[0]["PGA_roca_cm_s2"])
    target_rate = 1.0 / 475.0
    axis.axhline(
        target_rate,
        color="red",
        linestyle="--",
        linewidth=1.2,
        label=r"$\lambda=1/475$ año$^{-1}$",
    )
    axis.axvline(
        design_pga,
        color="red",
        linestyle="--",
        linewidth=1.2,
    )
    axis.annotate(
        "PGA (Tr = 475 años) = "
        f"{design_pga:.1f} cm/s² ({design_pga / 981.0:.3f} g)",
        xy=(design_pga, target_rate),
        xytext=(300.0, 0.025),
        color="red",
        fontsize=9,
        arrowprops={"arrowstyle": "->", "color": "red", "linewidth": 1.2},
        bbox={"boxstyle": "round,pad=0.3", "fc": "white", "ec": "red", "alpha": 0.9},
    )
    axis.set_yscale("log")
    axis.set_xlim(0.0, PGA_HAZARD_FIGURE_MAX_CM_S2)
    axis.set_ylim(1e-5, 1e1)
    axis.set_xlabel("PGA en roca (cm/s²)")
    axis.set_ylabel("Tasa anual de excedencia, λ (año⁻¹)")
    axis.set_title(
        "Amenaza sísmica probabilística en roca — Neiva\n"
        "Zhao et al. (2006), Clase I"
    )
    axis.grid(True, which="both", alpha=0.25)
    axis.legend(frameon=False, loc="upper right")
    figure.tight_layout()
    path = FIGURE_DIRECTORY / "curva_amenaza_neiva_roca_pga_semi_log.png"
    figure.savefig(path, dpi=250, bbox_inches="tight")
    plt.close(figure)
    return path


def run_analysis() -> dict[str, pd.DataFrame]:
    catalog, configuration, parameters = _read_inputs()
    catalog["distancia_epicentral_km"] = _epicentral_distance_km(
        catalog["latitud"].to_numpy(dtype=float),
        catalog["longitud"].to_numpy(dtype=float),
        float(configuration["site_latitude"]),
        float(configuration["site_longitude"]),
    )
    catalog["distancia_epicentral_km"] = catalog[
        "distancia_epicentral_km"
    ].clip(lower=MINIMUM_DISTANCE_KM)

    pga_visible = np.arange(
        0.0,
        PGA_DISPLAY_MAX_CM_S2 + PGA_STEP_CM_S2 / 2.0,
        PGA_STEP_CM_S2,
    )
    pga_extended = np.arange(
        PGA_DISPLAY_MAX_CM_S2 + PGA_EXTENDED_STEP_CM_S2,
        PGA_CALCULATION_MAX_CM_S2 + PGA_EXTENDED_STEP_CM_S2 / 2.0,
        PGA_EXTENDED_STEP_CM_S2,
    )
    pga = np.concatenate((pga_visible, pga_extended))
    source_curves: dict[str, np.ndarray] = {}
    geometry_rows = []
    for source in RECURRENCE:
        source_events = catalog[catalog["fuente"] == source].copy()
        source_curves[source], geometry = _source_hazard(
            source,
            source_events,
            pga,
            parameters,
            configuration,
        )
        geometry_rows.append({"Sismofuente": source, **geometry})

    total_curve = np.sum(np.stack(list(source_curves.values())), axis=0)
    for curve in [*source_curves.values(), total_curve]:
        if not np.isfinite(curve).all() or np.any(curve < 0.0):
            raise ArithmeticError("La curva de amenaza contiene tasas no válidas.")
        if np.any(np.diff(curve) > 1e-10):
            raise ArithmeticError("La tasa de excedencia debe ser no creciente.")
    if not np.isclose(total_curve[0], sum(
        row["tasa_truncada_Mw_ge_3_por_anio"] for row in geometry_rows
    )):
        raise ArithmeticError("La tasa total en PGA=0 no conserva las tasas fuente.")

    design = _design_accelerations(pga, total_curve)
    fit_rows = [
        {
            "Sismofuente": source,
            **_fit_exponential(
                pga[pga <= PGA_DISPLAY_MAX_CM_S2],
                curve[pga <= PGA_DISPLAY_MAX_CM_S2],
            ),
        }
        for source, curve in source_curves.items()
    ]
    fit_rows.append(
        {
            "Sismofuente": "Amenaza total en roca",
            **_fit_exponential(
                pga[pga <= PGA_DISPLAY_MAX_CM_S2],
                total_curve[pga <= PGA_DISPLAY_MAX_CM_S2],
            ),
        }
    )
    fits = pd.DataFrame(fit_rows)

    DATA_DIRECTORY.mkdir(parents=True, exist_ok=True)
    FIGURE_DIRECTORY.mkdir(parents=True, exist_ok=True)
    hazard = pd.DataFrame(
        {
            "PGA_roca_cm_s2": pga,
            **{
                f"lambda_{source}_por_anio": curve
                for source, curve in source_curves.items()
            },
            "lambda_total_por_anio": total_curve,
        }
    )
    hazard.to_csv(
        DATA_DIRECTORY / "curva_amenaza_neiva_roca_paso4.csv",
        index=False,
        encoding="utf-8-sig",
    )
    design.to_csv(
        DATA_DIRECTORY / "pga_diseno_neiva_roca_periodos_retorno_paso4.csv",
        index=False,
        encoding="utf-8-sig",
    )
    fits.to_csv(
        DATA_DIRECTORY / "ajustes_exponenciales_curva_amenaza_paso4.csv",
        index=False,
        encoding="utf-8-sig",
    )
    geometry = pd.DataFrame(geometry_rows)
    geometry.to_csv(
        DATA_DIRECTORY / "geometria_empirica_sismofuentes_neiva_paso4.csv",
        index=False,
        encoding="utf-8-sig",
    )
    figure_path = _plot_hazard(pga, source_curves, total_curve, design)
    metadata = {
        "site": configuration["site_name"],
        "site_class": "I",
        "site_term_C1": 0.0,
        "intensity_measure": "PGA (cm/s2)",
        "sigma_ln": SIGMA_LN,
        "magnitude_minimum": MINIMUM_MAGNITUDE,
        "maximum_magnitudes": {
            source: values["Mw_max"] for source, values in RECURRENCE.items()
        },
        "recurrence_parameters": RECURRENCE,
        "distance_model": (
            "Empirical discrete epicentral-distance distribution by source "
            "from the classified 260-km catalog; Haversine to Parque Santander."
        ),
        "depth_model": (
            "Source median catalog depth, clipped to 0-125 km; depth adjustment "
            "is active for h >= 15 km."
        ),
        "magnitude_distribution": (
            "Truncated Gutenberg-Richter incremental rates over Mw 3.0 to "
            "the configured source Mw_max; 0.025-Mw bins."
        ),
        "hazard_definition": (
            "Sum over source and magnitude-bin rate mass of the mean lognormal "
            "conditional exceedance probability over empirical source distances."
        ),
        "exponential_fit": (
            "lambda = A*exp(-B*PGA), fitted by linear regression of ln(lambda) "
            "over 0-100 cm/s2; R2 calculated in ln(lambda) space."
        ),
        "pga_range_cm_s2": {
            "exponential_fit": [0.0, PGA_DISPLAY_MAX_CM_S2],
            "hazard_figure": [0.0, PGA_HAZARD_FIGURE_MAX_CM_S2],
            "calculation_for_return_periods": [
                0.0,
                PGA_CALCULATION_MAX_CM_S2,
            ],
        },
        "return_period_convention": "Annual exceedance rate lambda = 1/Tr.",
        "note": (
            "Mw_max values and sigma are user-specified scenarios. "
            "The exponential equations are descriptive fits; design PGA is "
            "interpolated directly on the numerical total hazard curve."
        ),
    }
    (DATA_DIRECTORY / "parametros_integracion_amenaza_paso4.json").write_text(
        json.dumps(metadata, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )

    print("PGA de diseño en roca (interpolación en la curva total):")
    print(design.to_string(index=False))
    print("\nAjustes exponenciales lambda=A*exp(-B*PGA):")
    print(fits.to_string(index=False))
    print(f"\nCurva: {DATA_DIRECTORY / 'curva_amenaza_neiva_roca_paso4.csv'}")
    print(f"Figura: {figure_path}")
    return {
        "curva": hazard,
        "diseno": design,
        "ajustes": fits,
        "geometria": geometry,
    }


if __name__ == "__main__":
    run_analysis()
