from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
CATALOG = (
    ROOT
    / "recursos"
    / "Qgis"
    / "datos"
    / "Sismos_260_clasificados_20261004_192107.csv"
)
TABLES = ROOT / "datos" / "datos estadistica"
FIGURES = ROOT / "figuras"
SOURCES = [
    "Cortical",
    "Subducción interplaca",
    "Subducción intraplaca",
]
EDGES = np.array([3.0, 4.0, 5.0, 6.0, 7.0, 8.0])
MIN_SEGMENT_YEARS = 10
MIN_SEGMENT_EVENTS = 3
MIN_POST_R2 = 0.90
SITE_LATITUDE = 2.9262767
SITE_LONGITUDE = -75.2892111
RADII = [
    ("Regional (260 km)", "regional", 260.0),
]
TABLE_RECURRENCE = TABLES / "tabla_recurrencia_gr_por_tramo.csv"
TABLE_FITS = TABLES / "ajuste_gutenberg_richter_sismofuentes.csv"


def _r_squared(observed: np.ndarray, predicted: np.ndarray) -> float | None:
    total = float(np.sum((observed - observed.mean()) ** 2))
    if total == 0:
        return None
    residual = float(np.sum((observed - predicted) ** 2))
    return 1.0 - residual / total


def _completeness(
    annual: np.ndarray, years: np.ndarray, lower_magnitude: float
) -> tuple[int | None, float | None, str]:
    n_events = int(annual.sum())
    if lower_magnitude >= 5.0:
        return int(years[0]), None, "Asumida completa: Mw >= 5.0"
    if n_events == 0:
        return None, None, "Sin eventos"

    cumulative = np.cumsum(annual)
    best: tuple[float, int, np.ndarray, np.ndarray] | None = None
    for split in range(
        MIN_SEGMENT_YEARS, len(years) - MIN_SEGMENT_YEARS + 1
    ):
        pre_events = int(cumulative[split - 1])
        post_events = n_events - pre_events
        if min(pre_events, post_events) < MIN_SEGMENT_EVENTS:
            continue
        pre_fit = np.polyfit(years[:split], cumulative[:split], 1)
        post_fit = np.polyfit(years[split:], cumulative[split:], 1)
        pre_prediction = np.polyval(pre_fit, years[:split])
        post_prediction = np.polyval(post_fit, years[split:])
        error = float(
            np.sum((cumulative[:split] - pre_prediction) ** 2)
            + np.sum((cumulative[split:] - post_prediction) ** 2)
        )
        if best is None or error < best[0]:
            best = (error, split, pre_fit, post_fit)

    if best is None:
        return int(years[0]), None, "Referencia completa: quiebre no estimable"

    _, split, _, post_fit = best
    post_prediction = np.polyval(post_fit, years[split:])
    post_r2 = _r_squared(cumulative[split:], post_prediction)
    if post_r2 is None or post_r2 < MIN_POST_R2:
        return int(years[0]), post_r2, "Referencia completa: R2_post < 0.90"
    return int(years[split]), post_r2, "Estimada"


def _catalog() -> pd.DataFrame:
    data = pd.read_csv(CATALOG, encoding="latin1")
    data["Mw"] = pd.to_numeric(data["Mw"], errors="coerce")
    data["latitud"] = pd.to_numeric(data["latitud"], errors="coerce")
    data["longitud"] = pd.to_numeric(data["longitud"], errors="coerce")
    dates = pd.to_datetime(data["hora"], format="mixed", errors="coerce", utc=True)
    required = data[["Mw", "latitud", "longitud", "fuente"]]
    if required.isna().any().any() or dates.isna().any():
        raise ValueError("El catálogo contiene campos requeridos vacíos.")
    if not data["fuente"].isin(SOURCES).all():
        raise ValueError("El catálogo contiene una sismofuente no configurada.")
    data["year"] = dates.dt.year
    last_complete_year = int(dates.max().year)
    if (dates.max().month, dates.max().day) != (12, 31):
        last_complete_year -= 1
    data = data[data["year"] <= last_complete_year].copy()
    lat = np.radians(data["latitud"].to_numpy(dtype=float))
    lon = np.radians(data["longitud"].to_numpy(dtype=float))
    site_lat = np.radians(SITE_LATITUDE)
    site_lon = np.radians(SITE_LONGITUDE)
    haversine = (
        np.sin((lat - site_lat) / 2) ** 2
        + np.cos(site_lat)
        * np.cos(lat)
        * np.sin((lon - site_lon) / 2) ** 2
    )
    data["distance_km"] = (
        6371.0 * 2 * np.arcsin(np.sqrt(np.clip(haversine, 0.0, 1.0)))
    )
    return data


def _table_for_radius(
    catalog: pd.DataFrame, scale_name: str, scale_slug: str, radius: float
) -> pd.DataFrame:
    data = catalog[catalog["distance_km"] <= radius].copy()
    years = np.arange(int(catalog["year"].min()), int(catalog["year"].max()) + 1)
    rows = []
    for source in SOURCES:
        source_data = data[data["fuente"] == source]
        for lower, upper in zip(EDGES[:-1], EDGES[1:]):
            events = source_data[source_data["Mw"].between(
                lower, upper, inclusive="left"
            )]
            counts = (
                events.groupby("year")
                .size()
                .reindex(years, fill_value=0)
                .to_numpy(dtype=int)
            )
            start_year, post_r2, status = _completeness(
                counts, years, float(lower)
            )
            if start_year is None:
                duration = np.nan
                rate = 0.0
            else:
                duration = int(years[-1] - start_year + 1)
                complete_events = events[events["year"] >= start_year]
                rate = len(complete_events) / duration
            rows.append(
                {
                    "Escala": scale_name,
                    "Sismofuente": source,
                    "Intervalo Mw": f"[{lower:.0f}, {upper:.0f})",
                    "N": len(events),
                    "Año inicio completo": start_year,
                    "Delta t completo (años)": duration,
                    "R2 post": post_r2,
                    "lambda corregida (eventos/año)": rate,
                    "Estado": status,
                }
            )
    result = pd.DataFrame(rows)
    result = result[result["N"] > 0].copy()
    result.to_csv(
        TABLES / f"completez_agregada_fuentes_{scale_slug}.csv",
        index=False,
        encoding="utf-8-sig",
        float_format="%.6f",
    )
    return result


def _plot_paired(tables: dict[str, pd.DataFrame]) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.2), sharex=True)
    colors = ["#1f77b4", "#9467bd", "#e377c2"]
    intervals = [f"[{low:.0f}, {high:.0f})" for low, high in zip(EDGES[:-1], EDGES[1:])]
    x = np.arange(len(intervals))
    width = 0.24
    table = tables["regional"]
    ax_n, ax_rate = axes
    for index, (source, color) in enumerate(zip(SOURCES, colors)):
        subset = (
            table[table["Sismofuente"] == source]
            .set_index("Intervalo Mw")
            .reindex(intervals)
        )
        counts = subset["N"].to_numpy(dtype=float)
        rates = subset["lambda corregida (eventos/año)"].to_numpy(dtype=float)
        positions = x + (index - 1) * width
        ax_n.bar(
            positions,
            counts,
            width=width,
            color=color,
            label=source,
        )
        ax_rate.plot(
            x,
            rates,
            marker="o",
            linewidth=1.8,
            markersize=4,
            color=color,
            label=source,
        )
    ax_n.set_title("A. Eventos por intervalo — Regional (260 km)")
    ax_rate.set_title("B. Tasa anual corregida — Regional (260 km)")
    ax_n.set_yscale("log")
    ax_n.set_ylabel("Número de eventos, N")
    ax_n.set_xticks(x, intervals)
    ax_n.grid(axis="y", which="both", alpha=0.25)
    ax_n.legend(fontsize=8, loc="best")
    ax_rate.set_yscale("log")
    ax_rate.set_ylabel(r"$\lambda$ (eventos/año)")
    ax_rate.set_xticks(x, intervals)
    ax_rate.grid(axis="y", which="both", alpha=0.25)
    ax_rate.legend(fontsize=8, loc="best")
    fig = axes[0].figure
    fig.suptitle("Completez y recurrencia por sismofuente (ΔMw = 1.0)")
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(
        FIGURES / "completez_recurrencia_fuentes_intervalos_1mw.png",
        dpi=220,
        bbox_inches="tight",
    )
    plt.close(fig)


def _plot_regional_recurrence() -> None:
    recurrence = pd.read_csv(TABLE_RECURRENCE, encoding="utf-8-sig")
    fits = pd.read_csv(TABLE_FITS, encoding="utf-8-sig")
    scale = "Regional (260 km)"
    fig, ax = plt.subplots(figsize=(10.5, 6.5))
    colors = {
        "Cortical": "#1f77b4",
        "Subducción interplaca": "#9467bd",
        "Subducción intraplaca": "#e377c2",
    }
    x_fit = np.linspace(3.0, 7.0, 200)
    for source, color in colors.items():
        points = recurrence[
            recurrence["Escala"].eq(scale)
            & recurrence["Conjunto"].eq(source)
            & recurrence["log10(lambda)"].notna()
            & recurrence["Base tasa"].str.startswith("Completez", na=False)
        ]
        fit = fits[
            fits["Escala"].eq(scale) & fits["Conjunto"].eq(source)
        ]
        if points.empty or fit.empty:
            raise ValueError(f"Faltan datos de recurrencia regional para {source}.")
        fit_row = fit.iloc[0]
        a = float(fit_row["a"])
        b = float(fit_row["b"])
        r_squared = float(fit_row["R2"])
        ax.scatter(
            points["Mw mínimo"],
            points["log10(lambda)"],
            marker="s",
            color=color,
            label=f"{source}: tasas corregidas",
        )
        ax.plot(
            x_fit,
            a - b * x_fit,
            color=color,
            linewidth=2,
            label=(
                f"{source}: y = -{b:.3f}x + {a:.3f}; "
                f"R² = {r_squared:.3f}"
            ),
        )
    ax.set_title("Recurrencia Gutenberg–Richter — Regional (260 km)")
    ax.set_xlabel(r"Magnitud umbral $M_w$")
    ax.set_ylabel(r"$\log_{10}\lambda(M_w\geq M)$ (eventos/año)")
    ax.set_xlim(3.0, 7.0)
    ax.grid(which="both", alpha=0.25)
    ax.legend(fontsize=8, loc="best")
    fig.tight_layout()
    fig.savefig(
        FIGURES / "curvas_recurrencia_sismofuentes_regional.png",
        dpi=220,
        bbox_inches="tight",
    )
    plt.close(fig)


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    catalog = _catalog()
    tables = {
        slug: _table_for_radius(catalog, name, slug, radius)
        for name, slug, radius in RADII
    }
    _plot_paired(tables)
    _plot_regional_recurrence()
    for table in tables.values():
        print(table.to_string(index=False))


if __name__ == "__main__":
    main()
