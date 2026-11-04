from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D


DELIVERY_ROOT = Path(__file__).resolve().parents[1]
CATALOG_PATH = (
    DELIVERY_ROOT
    / "recursos"
    / "Qgis"
    / "datos"
    / "Sismos_260_clasificados_20261004_192107.csv"
)
OUTPUT_DIRECTORY = DELIVERY_ROOT / "datos" / "datos estadistica"
FIGURES_DIRECTORY = DELIVERY_ROOT / "figuras"

CONFIG: dict[str, Any] = {
    "site_name": "Neiva - Parque Santander",
    "site_latitude": 2.9262767,
    "site_longitude": -75.2892111,
    "site_class": "III",
    "intensity_measure": "PGA",
    "source_type_mapping": {
        "Cortical": "crustal",
        "Subducción interplaca": "interface",
        "Subducción intraplaca": "slab",
    },
    "reverse_fault_by_source": {
        "Cortical": False,
        "Subducción interplaca": False,
        "Subducción intraplaca": False,
    },
    "magnitude_edges": [
        3.0,
        3.3,
        3.5,
        4.0,
        4.5,
        5.0,
        5.5,
        6.0,
        6.5,
        7.0,
        7.5,
    ],
    "reference_magnitudes": [3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0],
    "reference_depth": "source_median",
    "depth_cap_km": 125.0,
    "depth_reference_km": 15.0,
    "minimum_distance_km": 0.1,
    "earth_radius_km": 6371.0,
}

MODEL_PARAMETERS: dict[str, dict[str, Any]] = {
    "PGA": {
        "label": "PGA",
        "period_s": 0.0,
        "a": 1.101,
        "b": -0.00564,
        "c": 0.0055,
        "d": 1.080,
        "e": 0.01412,
        "reverse_fault": 0.251,
        "interface": 0.000,
        "slab": 2.607,
        "slab_distance": -0.528,
        "site_terms": {
            "I": 1.111,
            "II": 1.344,
            "III": 1.355,
            "IV": 1.420,
        },
    },
    "Sa(0.5 s)": {
        "label": "Sa(0.5 s)",
        "period_s": 0.5,
        "a": 1.250,
        "b": -0.00338,
        "c": 0.0060,
        "d": 1.008,
        "e": 0.01114,
        "reverse_fault": 0.247,
        "interface": -0.053,
        "slab": 2.629,
        "slab_distance": -0.554,
        "site_terms": {
            "I": 0.071,
            "II": 0.515,
            "III": 0.934,
            "IV": 0.955,
        },
    },
}

def _validate_config(
    config: dict[str, Any], model_parameters: dict[str, dict[str, Any]]
) -> dict[str, Any]:
    required = {
        "site_name",
        "site_latitude",
        "site_longitude",
        "site_class",
        "intensity_measure",
        "source_type_mapping",
        "reverse_fault_by_source",
        "magnitude_edges",
        "reference_magnitudes",
        "reference_depth",
        "depth_cap_km",
        "depth_reference_km",
        "minimum_distance_km",
        "earth_radius_km",
    }
    missing = required.difference(config)
    if missing:
        raise ValueError(f"Faltan parámetros de configuración: {sorted(missing)}")
    if config["intensity_measure"] not in model_parameters:
        raise ValueError(
            "Medida de intensidad no configurada: "
            f"{config['intensity_measure']!r}"
        )
    parameters = model_parameters[config["intensity_measure"]]
    if config["site_class"] not in parameters["site_terms"]:
        raise ValueError(
            f"Clase de sitio no disponible: {config['site_class']!r}"
        )
    if not config["source_type_mapping"]:
        raise ValueError("Debe configurarse al menos una sismofuente.")
    unsupported_source_types = set(config["source_type_mapping"].values()) - {
        "crustal",
        "interface",
        "slab",
    }
    if unsupported_source_types:
        raise ValueError(
            f"Tipos tectónicos no soportados: {sorted(unsupported_source_types)}"
        )
    if any(
        not isinstance(enabled, bool)
        for enabled in config["reverse_fault_by_source"].values()
    ):
        raise ValueError("reverse_fault_by_source debe contener valores booleanos.")
    if not -90 <= float(config["site_latitude"]) <= 90:
        raise ValueError("La latitud del sitio debe estar entre -90 y 90 grados.")
    if not -180 <= float(config["site_longitude"]) <= 180:
        raise ValueError(
            "La longitud del sitio debe estar entre -180 y 180 grados."
        )
    edges = np.asarray(config["magnitude_edges"], dtype=float)
    if len(edges) < 2 or not np.all(np.diff(edges) > 0):
        raise ValueError("Los límites de magnitud deben ser estrictamente crecientes.")
    if float(config["minimum_distance_km"]) <= 0:
        raise ValueError("La distancia mínima debe ser positiva.")
    if float(config["earth_radius_km"]) <= 0:
        raise ValueError("El radio terrestre debe ser positivo.")
    if float(config["depth_cap_km"]) <= float(config["depth_reference_km"]):
        raise ValueError("La profundidad máxima debe superar la de referencia.")
    if config["reference_depth"] != "source_median":
        reference_depth = float(config["reference_depth"])
        if not 0 <= reference_depth <= float(config["depth_cap_km"]):
            raise ValueError(
                "La profundidad de referencia debe estar dentro de "
                "0 y depth_cap_km."
            )
    return parameters


def _epicentral_distance_km(
    latitude: np.ndarray,
    longitude: np.ndarray,
    site_latitude: float,
    site_longitude: float,
    earth_radius_km: float,
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
    return earth_radius_km * central_angle


def _source_terms(
    source: str,
    parameters: dict[str, Any],
    config: dict[str, Any],
) -> tuple[float, ...]:
    tectonic_type = config["source_type_mapping"].get(source)
    if tectonic_type is None:
        raise ValueError(f"No hay términos tectónicos configurados para {source!r}")
    reverse_fault = (
        float(parameters["reverse_fault"])
        if config["reverse_fault_by_source"].get(source, False)
        else 0.0
    )
    interface = 0.0
    slab = 0.0
    slab_distance = 0.0
    if tectonic_type == "interface":
        interface = float(parameters["interface"])
    elif tectonic_type == "slab":
        slab = float(parameters["slab"])
        slab_distance = float(parameters["slab_distance"])
    return reverse_fault, interface, slab, slab_distance


def _predict_intensity(
    magnitude: np.ndarray,
    distance_km: np.ndarray,
    depth_km: np.ndarray,
    source: str,
    parameters: dict[str, Any],
    config: dict[str, Any],
) -> np.ndarray:
    adjusted_depth = np.clip(
        depth_km, 0.0, float(config["depth_cap_km"])
    )
    distance_for_log = np.maximum(
        distance_km, float(config["minimum_distance_km"])
    )
    source_terms = _source_terms(source, parameters, config)
    depth_indicator = (
        adjusted_depth >= float(config["depth_reference_km"])
    ).astype(float)
    adjusted_distance = distance_for_log + float(parameters["c"]) * np.exp(
        float(parameters["d"]) * magnitude
    )
    ln_intensity = (
        float(parameters["a"]) * magnitude
        + float(parameters["b"]) * distance_for_log
        - np.log(adjusted_distance)
        + float(parameters["e"])
        * (adjusted_depth - float(config["depth_reference_km"]))
        * depth_indicator
        + source_terms[0]
        + source_terms[1]
        + source_terms[2]
        + source_terms[3] * np.log(distance_for_log)
        + float(parameters["site_terms"][config["site_class"]])
    )
    return np.exp(ln_intensity)


def _magnitude_label(lower: float, upper: float) -> str:
    return f"[{lower:g}, {upper:g})"


def _prepare_catalog(
    catalog_path: Path, config: dict[str, Any]
) -> pd.DataFrame:
    if not catalog_path.is_file():
        raise FileNotFoundError(f"No se encontró el catálogo: {catalog_path}")
    catalog = pd.read_csv(catalog_path, encoding="latin1")
    required_columns = {"latitud", "longitud", "Z", "Mw", "fuente"}
    missing_columns = required_columns.difference(catalog.columns)
    if missing_columns:
        raise ValueError(
            "Faltan columnas requeridas en el catálogo: "
            + ", ".join(sorted(missing_columns))
        )
    for column in ("latitud", "longitud", "Z", "Mw"):
        catalog[column] = pd.to_numeric(catalog[column], errors="coerce")
    if catalog[list(required_columns)].isna().any().any():
        raise ValueError(
            "El catálogo tiene valores vacíos o no numéricos en las variables "
            "necesarias para el modelo."
        )
    if not catalog["fuente"].isin(config["source_type_mapping"]).all():
        unexpected = sorted(
            set(catalog["fuente"]) - set(config["source_type_mapping"])
        )
        raise ValueError(f"Hay sismofuentes sin mapeo tectónico: {unexpected}")
    if ((catalog["latitud"] < -90) | (catalog["latitud"] > 90)).any():
        raise ValueError("El catálogo contiene latitudes fuera de rango.")
    if ((catalog["longitud"] < -180) | (catalog["longitud"] > 180)).any():
        raise ValueError("El catálogo contiene longitudes fuera de rango.")
    edges = np.asarray(config["magnitude_edges"], dtype=float)
    if (
        (catalog["Mw"] < edges[0]).any()
        or (catalog["Mw"] >= edges[-1]).any()
    ):
        raise ValueError(
            "Hay magnitudes fuera de magnitude_edges; amplíe los límites."
        )

    catalog["distancia_epicentral_km"] = _epicentral_distance_km(
        catalog["latitud"].to_numpy(dtype=float),
        catalog["longitud"].to_numpy(dtype=float),
        float(config["site_latitude"]),
        float(config["site_longitude"]),
        float(config["earth_radius_km"]),
    )
    catalog["profundidad_original_km"] = catalog["Z"]
    catalog["profundidad_ajustada_km"] = catalog["Z"].clip(
        lower=0.0, upper=float(config["depth_cap_km"])
    )
    catalog["profundidad_ajustada"] = (
        catalog["profundidad_original_km"]
        != catalog["profundidad_ajustada_km"]
    )
    distance_floor = float(config["minimum_distance_km"])
    catalog["distancia_usada_km"] = catalog[
        "distancia_epicentral_km"
    ].clip(lower=distance_floor)
    catalog["distancia_ajustada"] = (
        catalog["distancia_epicentral_km"] < distance_floor
    )

    magnitude_bins = pd.cut(
        catalog["Mw"],
        bins=edges,
        right=False,
        labels=[
            _magnitude_label(lower, upper)
            for lower, upper in zip(edges[:-1], edges[1:])
        ],
        include_lowest=True,
    )
    if magnitude_bins.isna().any():
        raise ValueError("No se pudieron asignar todos los eventos a un intervalo.")
    catalog["intervalo_Mw"] = magnitude_bins.astype(str)
    return catalog


def _summarize_by_source_and_bin(
    events: pd.DataFrame, parameters: dict[str, Any], config: dict[str, Any]
) -> pd.DataFrame:
    edges = np.asarray(config["magnitude_edges"], dtype=float)
    rows = []
    for source, source_events in events.groupby("fuente", sort=True):
        for lower, upper in zip(edges[:-1], edges[1:]):
            subset = source_events[
                source_events["Mw"].between(lower, upper, inclusive="left")
            ]
            if subset.empty:
                continue
            representative_magnitude = float(subset["Mw"].mean())
            representative_depth = float(
                subset["profundidad_ajustada_km"].mean()
            )
            mean_distance = float(subset["distancia_epicentral_km"].mean())
            representative_intensity = float(
                _predict_intensity(
                    np.array([representative_magnitude]),
                    np.array([mean_distance]),
                    np.array([representative_depth]),
                    source,
                    parameters,
                    config,
                )[0]
            )
            predictions = subset["intensidad_predicha_cm_s2"]
            rows.append(
                {
                    "Sitio": config["site_name"],
                    "Sismofuente": source,
                    "Intervalo Mw": _magnitude_label(lower, upper),
                    "N eventos": int(len(subset)),
                    "Mw media": representative_magnitude,
                    "Profundidad media ajustada (km)": representative_depth,
                    "Distancia epicentral media (km)": mean_distance,
                    "Distancia epicentral mediana (km)": float(
                        subset["distancia_epicentral_km"].median()
                    ),
                    "intensidad modelo en R media (cm/s2)": representative_intensity,
                    "intensidad mediana por evento (cm/s2)": float(
                        predictions.median()
                    ),
                    "intensidad media geométrica por evento (cm/s2)": float(
                        np.exp(np.log(predictions).mean())
                    ),
                    "intensidad P16 entre eventos (cm/s2)": float(
                        predictions.quantile(0.16)
                    ),
                    "intensidad P84 entre eventos (cm/s2)": float(
                        predictions.quantile(0.84)
                    ),
                }
            )
    return pd.DataFrame(rows)


def _plot_attenuation_curves(
    events: pd.DataFrame,
    parameters: dict[str, Any],
    config: dict[str, Any],
) -> list[Path]:
    sources = list(config["source_type_mapping"])
    figure_paths = []
    magnitude_edges = np.asarray(config["magnitude_edges"], dtype=float)
    magnitude_cmap = plt.get_cmap("viridis")
    magnitude_norm = plt.Normalize(
        float(magnitude_edges[0]), float(magnitude_edges[-1])
    )
    reference_magnitudes = [
        float(magnitude)
        for magnitude in config["reference_magnitudes"]
        if events["Mw"].min() <= float(magnitude) <= events["Mw"].max()
    ]
    source_slugs = {
        "Cortical": "cortical",
        "Subducción interplaca": "interplaca",
        "Subducción intraplaca": "intraplaca",
    }
    for source in sources:
        subset = events[events["fuente"] == source]
        reference_depth = (
            float(subset["profundidad_ajustada_km"].median())
            if config["reference_depth"] == "source_median"
            else float(config["reference_depth"])
        )
        distance_min = max(
            float(subset["distancia_epicentral_km"].min()),
            float(config["minimum_distance_km"]),
        )
        distance_max = float(subset["distancia_epicentral_km"].max())
        distance_grid = np.geomspace(distance_min, distance_max, 160)
        source_magnitude_min = float(subset["Mw"].min())
        source_magnitude_max = float(subset["Mw"].max())
        for x_scale, suffix in (("linear", ""), ("log", "_logx")):
            fig, axis = plt.subplots(figsize=(10, 7))
            scatter = axis.scatter(
                subset["distancia_epicentral_km"],
                subset["intensidad_predicha_cm_s2"],
                c=subset["Mw"],
                cmap=magnitude_cmap,
                norm=magnitude_norm,
                s=22,
                alpha=0.42,
                edgecolors="none",
            )
            for magnitude in config["reference_magnitudes"]:
                if not source_magnitude_min <= float(magnitude) <= source_magnitude_max:
                    continue
                intensity_curve = _predict_intensity(
                    np.full_like(distance_grid, float(magnitude)),
                    distance_grid,
                    np.full_like(distance_grid, reference_depth),
                    source,
                    parameters,
                    config,
                )
                axis.plot(
                    distance_grid,
                    intensity_curve,
                    color=magnitude_cmap(magnitude_norm(float(magnitude))),
                    linewidth=1.2,
                    label=f"Mw={float(magnitude):g}",
                )
            axis.set_xscale(x_scale)
            axis.set_yscale("log")
            axis.set_xlabel("Distancia epicentral al sitio (km)")
            axis.set_ylabel(
                f"{parameters['label']} predicha (cm/s², mediana condicional)"
            )
            axis.set_title(source)
            axis.grid(True, which="both", alpha=0.25)
            axis.set_xlim(distance_min, distance_max)
            colorbar = fig.colorbar(scatter, ax=axis, pad=0.02)
            colorbar.set_label("Magnitud Mw")
            handles = [
                Line2D(
                    [0],
                    [0],
                    color=magnitude_cmap(magnitude_norm(magnitude)),
                    linewidth=2,
                    label=f"Mw={magnitude:g}",
                )
                for magnitude in reference_magnitudes
                if source_magnitude_min <= magnitude <= source_magnitude_max
            ]
            if handles:
                axis.legend(
                    handles=handles,
                    title="Curvas para magnitud fija",
                    loc="upper right",
                    framealpha=0.9,
                    fontsize=9,
                    title_fontsize=9,
                )
            fig.suptitle(
                f"{parameters['label']} según Zhao et al. (2006)\n"
                f"{config['site_name']} — {source}; clase "
                f"{config['site_class']}",
                fontsize=15,
            )
            fig.tight_layout(rect=(0, 0, 1, 0.93))
            figure_path = FIGURES_DIRECTORY / (
                f"curvas_atenuacion_neiva_{measure_slug(parameters)}_"
                f"{source_slugs.get(source, 'fuente')}{suffix}.png"
            )
            fig.savefig(figure_path, dpi=300, bbox_inches="tight")
            plt.close(fig)
            figure_paths.append(figure_path)
    return figure_paths


def measure_slug(parameters: dict[str, Any]) -> str:
    return "pga" if parameters["period_s"] == 0 else "sa_0p5s"


def run_analysis(
    config: dict[str, Any] | None = None,
    model_parameters: dict[str, dict[str, Any]] | None = None,
    catalog_path: Path | None = None,
) -> dict[str, pd.DataFrame]:
    active_config = dict(CONFIG if config is None else config)
    active_parameters = (
        MODEL_PARAMETERS if model_parameters is None else model_parameters
    )
    parameters = _validate_config(active_config, active_parameters)
    active_catalog_path = CATALOG_PATH if catalog_path is None else catalog_path
    events = _prepare_catalog(active_catalog_path, active_config)
    events["intensidad_predicha_cm_s2"] = np.nan
    for source, indexes in events.groupby("fuente").groups.items():
        source_events = events.loc[indexes]
        events.loc[indexes, "intensidad_predicha_cm_s2"] = _predict_intensity(
            source_events["Mw"].to_numpy(dtype=float),
            source_events["distancia_epicentral_km"].to_numpy(dtype=float),
            source_events["profundidad_original_km"].to_numpy(dtype=float),
            source,
            parameters,
            active_config,
        )
    if not np.isfinite(events["intensidad_predicha_cm_s2"]).all():
        raise ValueError("El modelo produjo predicciones no finitas.")
    summary = _summarize_by_source_and_bin(events, parameters, active_config)
    if summary.empty:
        raise ValueError("No hay grupos fuente-magnitud con eventos para resumir.")

    OUTPUT_DIRECTORY.mkdir(parents=True, exist_ok=True)
    FIGURES_DIRECTORY.mkdir(parents=True, exist_ok=True)
    measure_slug = (
        "pga"
        if active_config["intensity_measure"] == "PGA"
        else "sa_0p5s"
    )
    events_path = OUTPUT_DIRECTORY / (
        f"predicciones_atenuacion_neiva_{measure_slug}_por_evento.csv"
    )
    summary_path = OUTPUT_DIRECTORY / (
        f"resumen_atenuacion_neiva_{measure_slug}_por_fuente_magnitud.csv"
    )
    metadata_path = OUTPUT_DIRECTORY / (
        f"parametros_atenuacion_neiva_{measure_slug}.json"
    )
    events.to_csv(events_path, index=False, encoding="utf-8-sig")
    summary.to_csv(summary_path, index=False, encoding="utf-8-sig")
    figure_paths = _plot_attenuation_curves(
        events, parameters, active_config
    )

    depth_adjustments = int(events["profundidad_ajustada"].sum())
    distance_adjustments = int(events["distancia_ajustada"].sum())
    metadata = {
        "configuration": active_config,
        "model_parameters": parameters,
        "depth_adjustments_by_source": (
            events.groupby("fuente")["profundidad_ajustada"]
            .sum()
            .astype(int)
            .to_dict()
        ),
        "distance_adjustments": distance_adjustments,
        "assumptions": [
            "Site class III is an assumed medium-soil scenario, not a measured site profile.",
            "Epicentral distance is calculated from each event to the configured site.",
            "Negative catalog depths are set to 0 km and depths above the configured cap are capped.",
            "Reverse-fault, intra-event and inter-event terms are zero unless explicitly changed in the inputs.",
            "Catalog Mw values are used as supplied; this is not a local Colombian ground-motion model validation.",
        ],
        "predicted_quantity": (
            "Conditional median ground motion; inter-event and intra-event "
            "residuals set to zero."
        ),
        "magnitude_source": "Catalog Mw column; no additional conversion applied.",
        "reference": (
            "Zhao et al. (2006), as formulated and tabulated in the "
            "Entrega 2 reference document."
        ),
    }
    metadata_path.write_text(
        json.dumps(metadata, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    print(f"Sitio: {active_config['site_name']}")
    print(
        f"Coordenadas: {active_config['site_latitude']}, "
        f"{active_config['site_longitude']}"
    )
    print(
        f"Medida: {parameters['label']}; clase de sitio "
        f"{active_config['site_class']}; eventos: {len(events):,}"
    )
    print(
        f"Profundidades limitadas a 0-{active_config['depth_cap_km']:g} km: "
        f"{depth_adjustments}; distancias ajustadas al mínimo: "
        f"{distance_adjustments}."
    )
    print("\nResumen por sismofuente e intervalo:")
    print(summary.to_string(index=False))
    print(f"\nPredicciones por evento: {events_path}")
    print(f"Resumen: {summary_path}")
    print("Figuras:")
    for figure_path in figure_paths:
        print(f"  {figure_path}")
    print(f"Parámetros reproducibles: {metadata_path}")
    return {"eventos": events, "resumen": summary}


if __name__ == "__main__":
    run_analysis()
