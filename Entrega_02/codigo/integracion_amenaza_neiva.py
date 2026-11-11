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
# Desviación estándar total de Zhao et al. (2006) para PGA (Tabla 5, sigma_T).
SIGMA_LN = 0.723
# Término de sitio de Zhao et al. (2006) para PGA en roca, Clase I (Tabla 5).
SITE_TERM_ROCK_C1 = 1.111
# Sensibilidad a la magnitud máxima de cada fuente (+/- 0.5).
SENSITIVITY_CASES: tuple[dict[str, Any], ...] = (
    {"caso": "Caso base"},
    {"caso": "Mw_max menor en 0.5", "delta_mmax": -0.5},
    {"caso": "Mw_max mayor en 0.5", "delta_mmax": 0.5},
)
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
    depths_km: np.ndarray,
    source: str,
    parameters: dict[str, Any],
    depth_reference_km: float,
    depth_cap_km: float,
    site_term: float,
) -> np.ndarray:
    distance = np.maximum(distances_km, MINIMUM_DISTANCE_KM)
    depth = np.clip(depths_km, 0.0, depth_cap_km)
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
    depth_term = np.where(
        depth >= depth_reference_km,
        float(parameters["e"]) * (depth - depth_reference_km),
        0.0,
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
        + site_term
    )
    return np.exp(ln_median)


def _source_hazard(
    source: str,
    source_events: pd.DataFrame,
    pga: np.ndarray,
    parameters: dict[str, Any],
    configuration: dict[str, Any],
    site_term: float = SITE_TERM_ROCK_C1,
    sigma_ln: float = SIGMA_LN,
    delta_mmax: float = 0.0,
    distance_column: str = "distancia_hipocentral_km",
) -> tuple[np.ndarray, dict[str, float]]:
    recurrence = RECURRENCE[source]
    maximum_magnitude = recurrence["Mw_max"] + delta_mmax
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

    distances = source_events[distance_column].to_numpy(dtype=float)
    # Cada evento del catálogo aporta su propia pareja (distancia, profundidad).
    depths = source_events["Z"].to_numpy(dtype=float)
    reference_depth = float(
        np.clip(np.median(depths), 0.0, configuration["depth_cap_km"])
    )
    source_curve = np.zeros_like(pga)
    positive_pga = pga > 0.0
    log_pga = np.log(pga[positive_pga])
    for magnitude, annual_rate in zip(magnitude_midpoints, annual_rate_in_bins):
        median = _attenuation_median(
            float(magnitude),
            distances,
            depths,
            source,
            parameters,
            float(configuration["depth_reference_km"]),
            float(configuration["depth_cap_km"]),
            site_term,
        )
        standardized = (
            log_pga[:, np.newaxis] - np.log(median)[np.newaxis, :]
        ) / sigma_ln
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


DISAGGREGATION_MAGNITUDE_BIN = 0.1
DISAGGREGATION_DISTANCE_EDGES_KM = (0, 25, 50, 75, 100, 150, 200, np.inf)
SOURCE_COLORS = {
    "Cortical": "#1f77b4",
    "Subducción interplaca": "#ff7f0e",
    "Subducción intraplaca": "#2ca02c",
}


def _weighted_median(values: np.ndarray, weights: np.ndarray) -> float:
    order = np.argsort(values)
    cumulative = np.cumsum(weights[order])
    return float(values[order][np.searchsorted(cumulative, 0.5 * cumulative[-1])])


def _disaggregate(
    catalog: pd.DataFrame,
    target_pga: float,
    parameters: dict[str, Any],
    configuration: dict[str, Any],
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, float]]:
    """Desagrega la tasa de excedencia de target_pga en fuente, Mw y R.

    Cada contribución es Delta_nu(M_k) * (1/N_i) * P(PGA > y* | M_k, R_j, h_j),
    con la misma configuración que la curva de amenaza del caso base.
    """
    rows = []
    for source, recurrence in RECURRENCE.items():
        events = catalog[catalog["fuente"] == source]
        distances = events["distancia_hipocentral_km"].to_numpy(dtype=float)
        depths = events["Z"].to_numpy(dtype=float)
        edges = np.append(
            np.arange(MINIMUM_MAGNITUDE, recurrence["Mw_max"], MAGNITUDE_STEP),
            recurrence["Mw_max"],
        )
        edges = np.unique(np.round(edges, 6))
        for lower, upper in zip(edges[:-1], edges[1:]):
            magnitude = (lower + upper) / 2.0
            rate = 10.0 ** (recurrence["a"] - recurrence["b"] * lower) - 10.0 ** (
                recurrence["a"] - recurrence["b"] * upper
            )
            median = _attenuation_median(
                float(magnitude),
                distances,
                depths,
                source,
                parameters,
                float(configuration["depth_reference_km"]),
                float(configuration["depth_cap_km"]),
                SITE_TERM_ROCK_C1,
            )
            probability = ndtr(-(np.log(target_pga) - np.log(median)) / SIGMA_LN)
            contribution = rate * probability / len(distances)
            rows.append(
                pd.DataFrame(
                    {
                        "Sismofuente": source,
                        "Mw": magnitude,
                        "R_hipocentral_km": distances,
                        "profundidad_km": np.clip(depths, 0.0, None),
                        "contribucion_anual": contribution,
                    }
                )
            )
    detail = pd.concat(rows, ignore_index=True)
    total = float(detail["contribucion_anual"].sum())
    detail["fraccion"] = detail["contribucion_anual"] / total

    summary_rows = []
    for source, group in [*detail.groupby("Sismofuente", sort=False), ("Total", detail)]:
        weights = group["fraccion"].to_numpy()
        share = float(weights.sum())
        normalized = weights / share
        summary_rows.append(
            {
                "Sismofuente": source,
                "Aporte_pct": 100.0 * share,
                "Mw_media": float(np.sum(normalized * group["Mw"])),
                "R_media_km": float(np.sum(normalized * group["R_hipocentral_km"])),
                "R_mediana_km": _weighted_median(
                    group["R_hipocentral_km"].to_numpy(), weights
                ),
                "profundidad_media_km": float(
                    np.sum(normalized * group["profundidad_km"])
                ),
            }
        )
    summary = pd.DataFrame(summary_rows)

    magnitude_bins = (
        np.floor(detail["Mw"] / DISAGGREGATION_MAGNITUDE_BIN + 1e-9)
        * DISAGGREGATION_MAGNITUDE_BIN
    )
    by_magnitude = detail.groupby(magnitude_bins)["fraccion"].sum()
    extra = {
        "lambda_total_por_anio": total,
        "Mw_modal": float(by_magnitude.idxmax() + DISAGGREGATION_MAGNITUDE_BIN / 2),
    }
    return detail, summary, extra


def _plot_disaggregation(detail: pd.DataFrame, target_pga: float) -> Path:
    figure, (axis_m, axis_r) = plt.subplots(1, 2, figsize=(12, 4.8), dpi=300)
    magnitude_bins = (
        np.floor(detail["Mw"] / DISAGGREGATION_MAGNITUDE_BIN + 1e-9)
        * DISAGGREGATION_MAGNITUDE_BIN
    ).round(2)
    distance_edges = np.asarray(DISAGGREGATION_DISTANCE_EDGES_KM, dtype=float)
    distance_labels = [
        f"≥{int(lower)}" if np.isinf(upper) else f"{int(lower)}-{int(upper)}"
        for lower, upper in zip(distance_edges[:-1], distance_edges[1:])
    ]
    distance_bins = pd.cut(
        detail["R_hipocentral_km"],
        bins=distance_edges,
        labels=distance_labels,
        right=False,
    )
    all_magnitudes = np.sort(magnitude_bins.unique())
    bottom_m = np.zeros(len(all_magnitudes))
    bottom_r = np.zeros(len(distance_labels))
    for source, color in SOURCE_COLORS.items():
        mask = detail["Sismofuente"] == source
        per_m = (
            detail.loc[mask, "fraccion"]
            .groupby(magnitude_bins[mask])
            .sum()
            .reindex(all_magnitudes, fill_value=0.0)
            .to_numpy()
            * 100.0
        )
        axis_m.bar(
            all_magnitudes + DISAGGREGATION_MAGNITUDE_BIN / 2,
            per_m,
            width=DISAGGREGATION_MAGNITUDE_BIN * 0.9,
            bottom=bottom_m,
            color=color,
            label=source,
        )
        bottom_m += per_m
        per_r = (
            detail.loc[mask, "fraccion"]
            .groupby(distance_bins[mask], observed=False)
            .sum()
            .reindex(distance_labels, fill_value=0.0)
            .to_numpy()
            * 100.0
        )
        axis_r.bar(distance_labels, per_r, bottom=bottom_r, color=color, label=source)
        bottom_r += per_r
    axis_m.set_xlabel("Magnitud Mw")
    axis_m.set_ylabel("Contribución a λ(475) (%)")
    axis_m.set_title("Desagregación en magnitud")
    axis_m.legend(frameon=False, loc="upper left")
    axis_r.set_xlabel("Distancia hipocentral (km)")
    axis_r.set_ylabel("Contribución a λ(475) (%)")
    axis_r.set_title("Desagregación en distancia")
    axis_m.set_ylim(0.0, 1.12 * bottom_m.max())
    axis_r.set_ylim(0.0, 1.12 * bottom_r.max())
    for axis in (axis_m, axis_r):
        axis.grid(True, axis="y", alpha=0.3)
    figure.suptitle(
        f"Desagregación para Tr = 475 años (PGA = {target_pga / 981.0:.3f} g, roca)"
    )
    figure.tight_layout()
    path = FIGURE_DIRECTORY / "desagregacion_neiva_475.png"
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
    # Distancia hipocentral: sqrt(Repi^2 + h^2), con h >= 0 (sin el tope de
    # 125 km, que solo aplica al término de profundidad de Zhao).
    catalog["distancia_hipocentral_km"] = np.hypot(
        catalog["distancia_epicentral_km"], catalog["Z"].clip(lower=0.0)
    ).clip(lower=MINIMUM_DISTANCE_KM)

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

    sensitivity_rows = []
    for case in SENSITIVITY_CASES:
        case_total = np.zeros_like(pga)
        for source in RECURRENCE:
            curve, _ = _source_hazard(
                source,
                catalog[catalog["fuente"] == source],
                pga,
                parameters,
                configuration,
                site_term=case.get("site_term", SITE_TERM_ROCK_C1),
                sigma_ln=case.get("sigma", SIGMA_LN),
                delta_mmax=case.get("delta_mmax", 0.0),
                distance_column=case.get(
                    "distance_column", "distancia_hipocentral_km"
                ),
            )
            case_total += curve
        target = 1.0 / 475.0
        if case_total[-1] > target:
            case_pga = float("nan")  # fuera del rango calculado
        else:
            case_pga = float(
                np.interp(np.log(target), np.log(case_total[::-1]), pga[::-1])
            )
        sensitivity_rows.append(
            {
                "Caso": case["caso"],
                "PGA_475_cm_s2": case_pga,
                "PGA_475_g": case_pga / 981.0,
            }
        )
    sensitivity = pd.DataFrame(sensitivity_rows)
    sensitivity.to_csv(
        DATA_DIRECTORY / "sensibilidad_pga_475_paso4.csv",
        index=False,
        encoding="utf-8-sig",
    )

    metadata = {
        "site": configuration["site_name"],
        "site_class": "I",
        "site_term_C1": SITE_TERM_ROCK_C1,
        "intensity_measure": "PGA (cm/s2)",
        "sigma_ln": SIGMA_LN,
        "magnitude_minimum": MINIMUM_MAGNITUDE,
        "maximum_magnitudes": {
            source: values["Mw_max"] for source, values in RECURRENCE.items()
        },
        "recurrence_parameters": RECURRENCE,
        "distance_model": (
            "Empirical discrete hypocentral-distance distribution by source "
            "from the classified 260-km catalog: sqrt(Repi^2 + h^2), with "
            "Repi by Haversine to Parque Santander."
        ),
        "depth_model": (
            "Each catalog event uses its own depth, clipped to 0-125 km for "
            "the Zhao depth term; depth adjustment is active for h >= 15 km."
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
    design_pga_475 = float(
        design.loc[design["Periodo_retorno_anios"] == 475, "PGA_roca_cm_s2"].iloc[0]
    )
    disaggregation_detail, disaggregation, disaggregation_extra = _disaggregate(
        catalog, design_pga_475, parameters, configuration
    )
    disaggregation.to_csv(
        DATA_DIRECTORY / "desagregacion_475_paso4.csv",
        index=False,
        encoding="utf-8-sig",
    )
    disaggregation_figure = _plot_disaggregation(
        disaggregation_detail, design_pga_475
    )
    print("\nDesagregación para Tr = 475 años:")
    print(disaggregation.to_string(index=False))
    print(disaggregation_extra)
    print(f"Figura: {disaggregation_figure}")

    print("\nSensibilidad del PGA para Tr = 475 años:")
    print(sensitivity.to_string(index=False))
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
