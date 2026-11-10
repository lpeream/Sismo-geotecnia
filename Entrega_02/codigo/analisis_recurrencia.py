from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import unicodedata

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


DELIVERY_ROOT = Path(__file__).resolve().parents[1]
CATALOG_PATH = (
    DELIVERY_ROOT
    / "recursos"
    / "Qgis"
    / "datos"
    / "Sismos_260_clasificados_20261004_192107.csv"
)
MAGNITUDE_EDGES = np.array(
    [3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5]
)
MIN_SEGMENT_YEARS = 10
MIN_SEGMENT_EVENTS = 3
MIN_POST_R2 = 0.90
MIN_FIT_BINS = 3
FULL_CATALOG_FROM_MAGNITUDE = 4.5
NEIVA_LATITUDE = 2.9262767
NEIVA_LONGITUDE = -75.2892111
LOCAL_RADIUS_KM = 50.0
REGIONAL_RADIUS_KM = 260.0


@dataclass
class CompletenessFit:
    start_year: int | None
    n_events: int
    breakpoint_r2: float | None
    pre_r2: float | None
    post_r2: float | None
    pre_slope: float | None
    post_slope: float | None
    status: str
    years: np.ndarray
    cumulative: np.ndarray
    split_index: int | None
    pre_prediction: np.ndarray | None
    post_prediction: np.ndarray | None


def _r_squared(observed: np.ndarray, predicted: np.ndarray) -> float | None:
    total = float(np.sum((observed - observed.mean()) ** 2))
    if total == 0:
        return None
    residual = float(np.sum((observed - predicted) ** 2))
    return 1.0 - residual / total


def _fit_completeness(
    events: pd.DataFrame, years: np.ndarray
) -> CompletenessFit:
    annual_counts = (
        events.groupby("year").size().reindex(years, fill_value=0).to_numpy()
    )
    cumulative = np.cumsum(annual_counts)
    n_events = int(annual_counts.sum())

    def unestimated(status: str) -> CompletenessFit:
        return CompletenessFit(
            None,
            n_events,
            None,
            None,
            None,
            None,
            None,
            status,
            years,
            cumulative,
            None,
            None,
            None,
        )

    if n_events == 0:
        return unestimated("Sin eventos")

    best: tuple[float, int, np.ndarray, np.ndarray] | None = None
    for split in range(
        MIN_SEGMENT_YEARS, len(years) - MIN_SEGMENT_YEARS + 1
    ):
        pre_events = int(cumulative[split - 1])
        post_events = n_events - pre_events
        if (
            pre_events < MIN_SEGMENT_EVENTS
            or post_events < MIN_SEGMENT_EVENTS
        ):
            continue

        pre_x = years[:split]
        post_x = years[split:]
        pre_fit = np.polyfit(pre_x, cumulative[:split], 1)
        post_fit = np.polyfit(post_x, cumulative[split:], 1)
        pre_prediction = np.polyval(pre_fit, pre_x)
        post_prediction = np.polyval(post_fit, post_x)
        error = float(
            np.sum((cumulative[:split] - pre_prediction) ** 2)
            + np.sum((cumulative[split:] - post_prediction) ** 2)
        )
        if best is None or error < best[0]:
            best = (error, split, pre_fit, post_fit)

    if best is None:
        return unestimated("Sin quiebre; usar tasa de referencia sin corregir")

    _, split, pre_fit, post_fit = best
    pre_prediction = np.polyval(pre_fit, years[:split])
    post_prediction = np.polyval(post_fit, years[split:])
    pre_r2 = _r_squared(cumulative[:split], pre_prediction)
    post_r2 = _r_squared(cumulative[split:], post_prediction)
    full_prediction = np.concatenate([pre_prediction, post_prediction])
    breakpoint_r2 = _r_squared(cumulative, full_prediction)
    if post_r2 is None or post_r2 < MIN_POST_R2:
        status = f"No validado: R² posterior < {MIN_POST_R2:.2f}"
        start_year = None
    elif n_events < 10:
        status = f"Aproximado: N={n_events}; muestra baja"
        start_year = int(years[split])
    else:
        status = "Estimado"
        start_year = int(years[split])

    return CompletenessFit(
        start_year,
        n_events,
        breakpoint_r2,
        pre_r2,
        post_r2,
        float(pre_fit[0]),
        float(post_fit[0]),
        status,
        years,
        cumulative,
        split,
        pre_prediction,
        post_prediction,
    )


def _magnitude_bins() -> list[tuple[float, float, str]]:
    return [
        (float(lower), float(upper), f"[{lower:g}, {upper:g})")
        for lower, upper in zip(
            MAGNITUDE_EDGES[:-1], MAGNITUDE_EDGES[1:]
        )
    ]


def _make_completeness_tables(
    data: pd.DataFrame,
    years: np.ndarray,
    group_names: list[str],
    catalog_name: str,
    scale_name: str,
) -> tuple[pd.DataFrame, dict[tuple[str, str], CompletenessFit]]:
    rows = []
    fits: dict[tuple[str, str], CompletenessFit] = {}
    groups = [(catalog_name, data)] + [
        (source, data[data["source"] == source]) for source in group_names
    ]

    for group_name, group in groups:
        for lower, upper, label in _magnitude_bins():
            band = group[group["Mw"].between(lower, upper, inclusive="left")]
            if lower >= FULL_CATALOG_FROM_MAGNITUDE:
                annual_counts = (
                    band.groupby("year")
                    .size()
                    .reindex(years, fill_value=0)
                    .to_numpy()
                )
                fit = CompletenessFit(
                    int(years[0]),
                    int(annual_counts.sum()),
                    None,
                    None,
                    None,
                    None,
                    None,
                    "Asumido completo: catálogo total (Mw >= 4.5)",
                    years,
                    np.cumsum(annual_counts),
                    None,
                    None,
                    None,
                )
            else:
                fit = _fit_completeness(band, years)
            fits[(group_name, label)] = fit
            start_year = fit.start_year
            rows.append(
                {
                    "Escala": scale_name,
                    "Conjunto": group_name,
                    "Intervalo Mw": label,
                    "N total catálogo": fit.n_events,
                    "Año de inicio completo": start_year,
                    "t completo (años)": (
                        int(years[-1] - start_year + 1)
                        if start_year is not None
                        else np.nan
                    ),
                    "Pendiente prequiebre (eventos/año)": fit.pre_slope,
                    "Pendiente postquiebre (eventos/año)": fit.post_slope,
                    "R2 prequiebre": fit.pre_r2,
                    "R2 postquiebre": fit.post_r2,
                    "R2 segmentado": fit.breakpoint_r2,
                    "Estado": fit.status,
                }
            )
    return pd.DataFrame(rows), fits


def _save_annual_cumulative_counts(
    data: pd.DataFrame,
    years: np.ndarray,
    groups: list[str],
    output_path: Path,
    catalog_name: str,
    scale_name: str,
) -> None:
    rows = []
    all_groups = [(catalog_name, data)] + [
        (source, data[data["source"] == source]) for source in groups
    ]
    for group_name, group in all_groups:
        for lower, upper, label in _magnitude_bins():
            band = group[group["Mw"].between(lower, upper, inclusive="left")]
            counts = (
                band.groupby("year")
                .size()
                .reindex(years, fill_value=0)
                .to_numpy()
            )
            cumulative = np.cumsum(counts)
            rows.extend(
                {
                    "Escala": scale_name,
                    "Conjunto": group_name,
                    "Intervalo Mw": label,
                    "Año": int(year),
                    "N anual": int(count),
                    "N acumulado": int(cumulative_count),
                }
                for year, count, cumulative_count in zip(
                    years, counts, cumulative
                )
            )
    pd.DataFrame(rows).to_csv(
        output_path, index=False, encoding="utf-8-sig"
    )


def _save_completeness_figures(
    years: np.ndarray,
    groups: list[str],
    fits: dict[tuple[str, str], CompletenessFit],
    figures_dir: Path,
    catalog_name: str,
    scale_slug: str,
) -> None:
    all_groups = [catalog_name, *groups]
    for group_name in all_groups:
        labels = [
            label
            for _, _, label in _magnitude_bins()
            if fits[(group_name, label)].n_events > 0
        ]
        if not labels:
            continue
        n_rows = (len(labels) + 1) // 2
        fig, axes = plt.subplots(
            n_rows,
            2,
            figsize=(12, 3.2 * n_rows),
            sharex=True,
            squeeze=False,
        )
        for ax, label in zip(axes.flat, labels):
            fit = fits[(group_name, label)]
            ax.plot(
                years,
                fit.cumulative,
                color="#1f5d75",
                linewidth=1.5,
                label="Acumulado observado",
            )
            if fit.split_index is not None:
                split = fit.split_index
                ax.plot(
                    years[:split],
                    fit.pre_prediction,
                    color="#d8892c",
                    linestyle="--",
                    label="Ajuste prequiebre",
                )
                ax.plot(
                    years[split:],
                    fit.post_prediction,
                    color="#b24732",
                    linestyle="--",
                    label="Ajuste postquiebre",
                )
                ax.axvline(
                    years[split],
                    color="black",
                    linestyle=":",
                    linewidth=1,
                )
                if fit.start_year is not None:
                    subtitle = (
                        f"inicio={fit.start_year}, "
                        f"R² post={fit.post_r2:.3f}"
                    )
                else:
                    subtitle = (
                        f"sin validación; R² post="
                        f"{fit.post_r2:.3f}"
                    )
            else:
                subtitle = fit.status

            ax.set_title(f"{label} Mw; N={fit.n_events}; {subtitle}")
            ax.grid(alpha=0.25)
            ax.set_ylabel("N acumulado")
        for ax in axes.flat[len(labels) :]:
            ax.set_visible(False)
        last_row = (len(labels) - 1) // 2
        axes[last_row, 0].set_xlabel("Año")
        if len(labels) % 2 == 0:
            axes[last_row, 1].set_xlabel("Año")
        axes[0, 0].legend(fontsize=7, loc="upper left")
        fig.suptitle(
            f"Completez por intervalos de magnitud — {group_name}",
            y=0.995,
        )
        fig.tight_layout(rect=(0, 0, 1, 0.99))
        slug = _slug(group_name)
        if scale_slug == "local":
            slug = f"local_{slug}"
        fig.savefig(
            figures_dir / f"completez_acumulada_{slug}.png",
            dpi=220,
        )
        plt.close(fig)


def _make_recurrence_tables(
    data: pd.DataFrame,
    group_names: list[str],
    completeness: pd.DataFrame,
    first_complete_year: int,
    last_complete_year: int,
    catalog_name: str,
    scale_name: str,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    recurrence_rows = []
    fit_rows = []
    increment_rows = []
    groups = [(catalog_name, data)] + [
        (source, data[data["source"] == source]) for source in group_names
    ]

    for group_name, group in groups:
        valid_points = []
        for lower, upper, label in _magnitude_bins():
            completeness_row = completeness[
                (completeness["Conjunto"] == group_name)
                & (completeness["Intervalo Mw"] == label)
            ].iloc[0]
            start_year_value = completeness_row["Año de inicio completo"]
            status = completeness_row["Estado"]

            row: dict[str, object] = {
                "Escala": scale_name,
                "Conjunto": group_name,
                "Intervalo Mw": label,
                "Mw mínimo": lower,
                "Mw máximo": upper,
                "N incremental": np.nan,
                "N acumulado (misma ventana)": np.nan,
                "Año inicio completo": np.nan,
                "t (años)": np.nan,
                "lambda acumulada (año^-1)": np.nan,
                "log10(lambda)": np.nan,
                "Estado completez": status,
                "Base tasa": "",
            }
            if pd.notna(start_year_value):
                start_year = int(start_year_value)
                if str(status).startswith("Aproximado:"):
                    rate_basis = "Completez aproximada (N<10)"
                elif str(status).startswith("Asumido completo:"):
                    rate_basis = "Completez asumida: catálogo total"
                else:
                    rate_basis = "Completez estimada"
            else:
                start_year = first_complete_year
                rate_basis = "Referencia cruda; completez no estimable"
            duration = last_complete_year - start_year + 1
            observed = group[group["year"].between(start_year, last_complete_year)]
            if duration > 0:
                incremental = int(
                    observed["Mw"].between(
                        lower, upper, inclusive="left"
                    ).sum()
                )
                cumulative = int((observed["Mw"] >= lower).sum())
                # The cumulative count uses one consistent observation window:
                # the completeness interval determined for its lower threshold.
                crosscheck = sum(
                    int(
                        observed["Mw"].between(
                            band_lower, band_upper, inclusive="left"
                        ).sum()
                    )
                    for band_lower, band_upper, _ in _magnitude_bins()
                    if band_lower >= lower
                )
                if cumulative != crosscheck:
                    raise ValueError(
                        f"Los conteos incrementales no suman N acumulado "
                        f"para {group_name}, {label}"
                    )
                for increment_lower, increment_upper, increment_label in (
                    _magnitude_bins()
                ):
                    if increment_lower < lower:
                        continue
                    increment_rows.append(
                        {
                            "Escala": scale_name,
                            "Conjunto": group_name,
                            "Umbral acumulado Mw": label,
                            "Intervalo incremental Mw": increment_label,
                            "Año inicio ventana": start_year,
                            "t ventana (años)": duration,
                            "N incremental en ventana": int(
                                observed["Mw"].between(
                                    increment_lower,
                                    increment_upper,
                                    inclusive="left",
                                ).sum()
                            ),
                        }
                    )
                rate = cumulative / duration
                row.update(
                    {
                        "N incremental": incremental,
                        "N acumulado (misma ventana)": cumulative,
                        "Año inicio completo": start_year,
                        "t (años)": duration,
                        "lambda acumulada (año^-1)": rate,
                        "Base tasa": rate_basis,
                        "log10(lambda)": (
                            float(np.log10(rate)) if rate > 0 else np.nan
                        ),
                    }
                )
                if rate > 0 and pd.notna(start_year_value):
                    valid_points.append(
                        {"Mw mínimo": lower, "log10(lambda)": float(np.log10(rate))}
                    )
            else:
                row["Base tasa"] = "Sin exposición"
            recurrence_rows.append(row)

        if len(valid_points) >= MIN_FIT_BINS:
            x = np.array([point["Mw mínimo"] for point in valid_points])
            y = np.array([point["log10(lambda)"] for point in valid_points])
            slope, intercept = np.polyfit(x, y, 1)
            predicted = intercept + slope * x
            r_squared = _r_squared(y, predicted)
            fit = {
                "Escala": scale_name,
                "Conjunto": group_name,
                "a": float(intercept),
                "b": float(-slope),
                "R2": r_squared,
                "N intervalos ajuste": len(valid_points),
            }
            fit_rows.append(fit)
        else:
            fit = {
                "Escala": scale_name,
                "Conjunto": group_name,
                "a": np.nan,
                "b": np.nan,
                "R2": np.nan,
                "N intervalos ajuste": len(valid_points),
            }
            fit_rows.append(fit)

    recurrence = pd.DataFrame(recurrence_rows)
    fits = pd.DataFrame(fit_rows)
    increments = pd.DataFrame(increment_rows)
    for fit in fit_rows:
        selected = (
            recurrence["Escala"].eq(fit["Escala"])
            & recurrence["Conjunto"].eq(fit["Conjunto"])
        )
        recurrence.loc[selected, "a"] = fit["a"]
        recurrence.loc[selected, "b"] = fit["b"]
        recurrence.loc[selected, "R2 ajuste"] = fit["R2"]
    return recurrence, fits, increments


def _slug(value: str) -> str:
    normalized = unicodedata.normalize("NFKD", value)
    ascii_value = normalized.encode("ascii", "ignore").decode("ascii")
    return "_".join(ascii_value.lower().split())


def _make_uncorrected_recurrence(
    data: pd.DataFrame,
    groups: list[str],
    first_year: int,
    last_complete_year: int,
    catalog_name: str,
    scale_name: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    duration = last_complete_year - first_year + 1
    raw_rows = []
    fit_rows = []
    all_groups = [(catalog_name, data)] + [
        (source, data[data["source"] == source]) for source in groups
    ]
    for group_name, group in all_groups:
        valid = []
        for lower, _, label in _magnitude_bins():
            cumulative = int((group["Mw"] >= lower).sum())
            rate = cumulative / duration
            log_rate = float(np.log10(rate)) if rate > 0 else np.nan
            raw_rows.append(
                {
                    "Escala": scale_name,
                    "Conjunto": group_name,
                    "Intervalo Mw": label,
                    "Mw mínimo": lower,
                    "N acumulado sin corrección": cumulative,
                    "t ventana total (años)": duration,
                    "lambda sin corrección (año^-1)": rate,
                    "log10(lambda sin corrección)": log_rate,
                }
            )
            if rate > 0:
                valid.append((lower, log_rate))
        if len(valid) >= MIN_FIT_BINS:
            x = np.asarray([point[0] for point in valid])
            y = np.asarray([point[1] for point in valid])
            slope, intercept = np.polyfit(x, y, 1)
            fit_rows.append(
                {
                    "Escala": scale_name,
                    "Conjunto": group_name,
                    "a": float(intercept),
                    "b": float(-slope),
                    "R2": _r_squared(y, intercept + slope * x),
                    "N intervalos ajuste": len(valid),
                }
            )
        else:
            fit_rows.append(
                {
                    "Escala": scale_name,
                    "Conjunto": group_name,
                    "a": np.nan,
                    "b": np.nan,
                    "R2": np.nan,
                    "N intervalos ajuste": len(valid),
                }
            )
    return pd.DataFrame(raw_rows), pd.DataFrame(fit_rows)


def _save_source_comparison_figures(
    recurrence: pd.DataFrame,
    corrected_fits: pd.DataFrame,
    raw_recurrence: pd.DataFrame,
    raw_fits: pd.DataFrame,
    source_names: list[str],
    figures_dir: Path,
) -> None:
    scales = [("Local (50 km)", "local"), ("Regional (260 km)", "regional")]
    for scale_name, scale_slug in scales:
        for source in source_names:
            raw = raw_recurrence[
                raw_recurrence["Escala"].eq(scale_name)
                & raw_recurrence["Conjunto"].eq(source)
                & raw_recurrence["log10(lambda sin corrección)"].notna()
            ]
            corrected = recurrence[
                recurrence["Escala"].eq(scale_name)
                & recurrence["Conjunto"].eq(source)
                & recurrence["Base tasa"].str.startswith(
                    "Completez", na=False
                )
                & recurrence["log10(lambda)"].notna()
            ]
            reference = recurrence[
                recurrence["Escala"].eq(scale_name)
                & recurrence["Conjunto"].eq(source)
                & recurrence["Base tasa"].str.startswith(
                    "Referencia cruda", na=False
                )
                & recurrence["log10(lambda)"].notna()
            ]
            approximate = recurrence[
                recurrence["Escala"].eq(scale_name)
                & recurrence["Conjunto"].eq(source)
                & recurrence["Base tasa"].eq(
                    "Completez aproximada (N<10)"
                )
                & recurrence["log10(lambda)"].notna()
            ]
            if raw.empty and corrected.empty and approximate.empty:
                continue
            fig, ax = plt.subplots(figsize=(9, 6))
            if not raw.empty:
                ax.scatter(
                    raw["Mw mínimo"],
                    raw["log10(lambda sin corrección)"],
                    marker="o",
                    color="#d8892c",
                    label="Tasas sin corrección",
                    zorder=3,
                )
            if not corrected.empty:
                ax.scatter(
                    corrected["Mw mínimo"],
                    corrected["log10(lambda)"],
                    marker="s",
                    color="#1f5d75",
                    label="Tasas corregidas por completez",
                    zorder=4,
                )
            if not approximate.empty:
                ax.scatter(
                    approximate["Mw mínimo"],
                    approximate["log10(lambda)"],
                    marker="D",
                    color="#b24732",
                    label="Completez estimada con N<10; baja confianza",
                    zorder=5,
                )
            if not reference.empty:
                ax.scatter(
                    reference["Mw mínimo"],
                    reference["log10(lambda)"],
                    marker="x",
                    color="#777777",
                    label="Referencia cruda; completez no estimable",
                    zorder=2,
                )

            for fit_data, fit_kind, color, linestyle in (
                (raw_fits, "sin corrección", "#d8892c", "--"),
                (corrected_fits, "corregida", "#1f5d75", "-"),
            ):
                fit_rows = fit_data[
                    fit_data["Escala"].eq(scale_name)
                    & fit_data["Conjunto"].eq(source)
                ]
                if fit_rows.empty:
                    continue
                fit = fit_rows.iloc[0]
                if not np.isfinite(fit["a"]) or not np.isfinite(fit["b"]):
                    continue
                x_fit = np.linspace(
                    float(MAGNITUDE_EDGES[0]),
                    float(MAGNITUDE_EDGES[-2]),
                    200,
                )
                equation = (
                    f"{fit_kind}: y = -{fit['b']:.3f}*x + "
                    f"{fit['a']:.3f}; R² = {fit['R2']:.3f}"
                )
                ax.plot(
                    x_fit,
                    fit["a"] - fit["b"] * x_fit,
                    color=color,
                    linestyle=linestyle,
                    linewidth=1.8,
                    label=equation,
                )
            if raw.empty and corrected.empty:
                ax.text(
                    0.5,
                    0.5,
                    "Sin eventos en esta escala",
                    ha="center",
                    va="center",
                    transform=ax.transAxes,
                )
            elif corrected.empty:
                ax.text(
                    0.02,
                    0.03,
                    "Sin ajuste corregido: menos de tres intervalos "
                    "con completez estimada o asumida.",
                    ha="left",
                    va="bottom",
                    transform=ax.transAxes,
                    fontsize=8,
                )
            ax.set_xlabel(r"Magnitud umbral $M_w$")
            ax.set_ylabel(r"$y=\log_{10}\lambda(M_w\geq M)$ (eventos/año)")
            ax.set_title(f"{source} — {scale_name}")
            ax.grid(True, alpha=0.25)
            handles, labels = ax.get_legend_handles_labels()
            if handles:
                ax.legend(handles, labels, fontsize=8, loc="best")
            fig.tight_layout()
            fig.savefig(
                figures_dir
                / f"comparacion_recurrencia_cruda_completez_"
                f"{scale_slug}_{_slug(source)}.png",
                dpi=220,
            )
            plt.close(fig)


def _save_recurrence_plot(
    recurrence: pd.DataFrame,
    fits: pd.DataFrame,
    raw_recurrence: pd.DataFrame,
    figures_dir: Path,
) -> None:
    colors = {
        "Catálogo regional": "#303030",
        "Catálogo local": "#303030",
        "Cortical": "#1f77b4",
        "Subducción interplaca": "#9467bd",
        "Subducción intraplaca": "#e377c2",
    }
    max_threshold = float(MAGNITUDE_EDGES[-2])
    for scale_name, scale_slug in (
        ("Local (50 km)", "local"),
        ("Regional (260 km)", "regional"),
    ):
        scale_fits = fits[fits["Escala"].eq(scale_name)]
        positive_logs = pd.concat(
            [
                recurrence.loc[
                    recurrence["Escala"].eq(scale_name), "log10(lambda)"
                ].dropna(),
                raw_recurrence.loc[
                    raw_recurrence["Escala"].eq(scale_name),
                    "log10(lambda sin corrección)",
                ].dropna(),
            ],
            ignore_index=True,
        )
        finite_fits = scale_fits[
            scale_fits["Conjunto"].isin(colors)
            & scale_fits["a"].notna()
            & scale_fits["b"].notna()
        ]
        fit_endpoints = pd.concat(
            [
                finite_fits["a"] - finite_fits["b"] * MAGNITUDE_EDGES[0],
                finite_fits["a"] - finite_fits["b"] * max_threshold,
            ],
            ignore_index=True,
        )
        lowest_fit = (
            float(fit_endpoints.min()) if not fit_endpoints.empty else -1
        )
        lowest_data = (
            float(positive_logs.min()) if not positive_logs.empty else -1
        )
        plot_bottom = min(lowest_data - 0.35, lowest_fit - 0.2)
        zero_y = plot_bottom + 0.08
        fig, ax = plt.subplots(figsize=(10, 8))
        raw_legend_added = False
        corrected_legend_added = False
        has_zero_rates = False

        for fit in scale_fits.itertuples(index=False):
            color = colors.get(fit.Conjunto)
            if color is None:
                continue
            corrected = recurrence[
                recurrence["Escala"].eq(scale_name)
                & recurrence["Conjunto"].eq(fit.Conjunto)
                & recurrence["Base tasa"].str.startswith(
                    "Completez", na=False
                )
            ]
            raw = raw_recurrence[
                raw_recurrence["Escala"].eq(scale_name)
                & raw_recurrence["Conjunto"].eq(fit.Conjunto)
            ]
            raw_positive = raw[
                raw["log10(lambda sin corrección)"].notna()
            ]
            corrected_positive = corrected[
                corrected["log10(lambda)"].notna()
            ]
            if not raw_positive.empty:
                ax.scatter(
                    raw_positive["Mw mínimo"],
                    raw_positive["log10(lambda sin corrección)"],
                    marker="o",
                    facecolors="none",
                    edgecolors=color,
                    s=32,
                    label=(
                        "Sin corrección (círculos)"
                        if not raw_legend_added
                        else "_nolegend_"
                    ),
                    zorder=3,
                )
                raw_legend_added = True
            if not corrected_positive.empty:
                ax.scatter(
                    corrected_positive["Mw mínimo"],
                    corrected_positive["log10(lambda)"],
                    marker="s",
                    color=color,
                    s=27,
                    label=(
                        "Corregida (cuadrados)"
                        if not corrected_legend_added
                        else "_nolegend_"
                    ),
                    zorder=4,
                )
                corrected_legend_added = True

            raw_zeros = raw[raw["lambda sin corrección (año^-1)"].eq(0)]
            corrected_zeros = corrected[
                corrected["lambda acumulada (año^-1)"].eq(0)
            ]
            if not raw_zeros.empty:
                has_zero_rates = True
                ax.scatter(
                    raw_zeros["Mw mínimo"],
                    np.full(len(raw_zeros), zero_y),
                    marker="o",
                    facecolors="none",
                    edgecolors=color,
                    s=32,
                    zorder=3,
                )
            if not corrected_zeros.empty:
                has_zero_rates = True
                ax.scatter(
                    corrected_zeros["Mw mínimo"],
                    np.full(len(corrected_zeros), zero_y),
                    marker="s",
                    color=color,
                    s=27,
                    zorder=4,
                )

            if not np.isfinite(fit.a) or not np.isfinite(fit.b):
                continue
            x_fit = np.linspace(
                float(MAGNITUDE_EDGES[0]), max_threshold, 200
            )
            ax.plot(
                x_fit,
                fit.a - fit.b * x_fit,
                color=color,
                linewidth=1.8,
                label=(
                    f"{fit.Conjunto}: y = -{fit.b:.3f}*x + "
                    f"{fit.a:.3f}; R²={fit.R2:.3f}"
                ),
            )

        ax.set_xlabel(r"Magnitud umbral $M_w$", fontsize=13)
        ax.set_ylabel(
            r"$\log_{10}\lambda(M_w\geq M)$ (eventos/año)", fontsize=13
        )
        ax.tick_params(axis="both", labelsize=11)
        ax.set_xlim(float(MAGNITUDE_EDGES[0]), max_threshold)
        ax.set_ylim(bottom=plot_bottom)
        ax.set_title(
            f"Recurrencia Gutenberg--Richter — {scale_name}", fontsize=15
        )
        if has_zero_rates:
            fig.text(
                0.12,
                0.015,
                "Marcadores en el límite inferior: λ=0; log10(λ) no definido.",
                fontsize=9,
            )
        ax.grid(True, alpha=0.25)
        ax.legend(fontsize=10, loc="best")
        fig.tight_layout(rect=(0, 0.04, 1, 1))
        fig.savefig(
            figures_dir / f"curvas_recurrencia_sismofuentes_{scale_slug}.png",
            dpi=250,
        )
        plt.close(fig)


def run_analysis() -> dict[str, pd.DataFrame]:
    if not CATALOG_PATH.is_file():
        raise FileNotFoundError(f"No se encontró el catálogo: {CATALOG_PATH}")

    figures_dir = DELIVERY_ROOT / "figuras"
    tables_dir = DELIVERY_ROOT / "datos" / "datos estadistica"
    figures_dir.mkdir(parents=True, exist_ok=True)
    tables_dir.mkdir(parents=True, exist_ok=True)

    catalog = pd.read_csv(CATALOG_PATH, encoding="latin1")
    required = {"Mw", "hora", "fuente", "latitud", "longitud"}
    missing_columns = required.difference(catalog.columns)
    if missing_columns:
        raise ValueError(
            "Faltan columnas requeridas: " + ", ".join(sorted(missing_columns))
        )
    catalog["Mw"] = pd.to_numeric(catalog["Mw"], errors="coerce")
    catalog["latitud"] = pd.to_numeric(catalog["latitud"], errors="coerce")
    catalog["longitud"] = pd.to_numeric(catalog["longitud"], errors="coerce")
    catalog["date_utc"] = pd.to_datetime(
        catalog["hora"], format="mixed", errors="coerce", utc=True
    )
    if catalog[
        ["Mw", "date_utc", "fuente", "latitud", "longitud"]
    ].isna().any().any():
        raise ValueError(
            "El catálogo contiene magnitudes, fechas, fuentes o coordenadas vacías."
        )
    if (catalog["Mw"] < MAGNITUDE_EDGES[0]).any():
        raise ValueError(
            "El catálogo contiene magnitudes inferiores al primer límite "
            f"({MAGNITUDE_EDGES[0]:g}); revise los intervalos."
        )
    if (catalog["Mw"] >= MAGNITUDE_EDGES[-1]).any():
        raise ValueError(
            "El catálogo contiene magnitudes fuera del último límite "
            f"({MAGNITUDE_EDGES[-1]:g}); extienda MAGNITUDE_EDGES."
        )

    catalog["year"] = catalog["date_utc"].dt.year
    catalog["source"] = catalog["fuente"]
    latitude = np.radians(catalog["latitud"].to_numpy(dtype=float))
    longitude = np.radians(catalog["longitud"].to_numpy(dtype=float))
    site_latitude = np.radians(NEIVA_LATITUDE)
    site_longitude = np.radians(NEIVA_LONGITUDE)
    haversine = (
        np.sin((latitude - site_latitude) / 2) ** 2
        + np.cos(site_latitude)
        * np.cos(latitude)
        * np.sin((longitude - site_longitude) / 2) ** 2
    )
    catalog["distancia_neiva_km"] = (
        6371.0 * 2 * np.arcsin(np.sqrt(np.clip(haversine, 0.0, 1.0)))
    )
    last_timestamp = catalog["date_utc"].max()
    last_complete_year = int(last_timestamp.year)
    if (last_timestamp.month, last_timestamp.day) != (12, 31):
        last_complete_year -= 1
    data = catalog[catalog["year"] <= last_complete_year].copy()
    first_year = int(data["year"].min())
    years = np.arange(first_year, last_complete_year + 1)
    sources = sorted(data["fuente"].unique())

    print(f"Catálogo: {CATALOG_PATH}")
    print(f"Eventos totales: {len(catalog):,}")
    print(
        f"Eventos en años completos ({first_year}-{last_complete_year}): "
        f"{len(data):,}"
    )
    print(f"Fuentes: {', '.join(sources)}")
    print(f"Intervalos de magnitud: {', '.join(b[2] for b in _magnitude_bins())}")
    print(
        "Quiebre: mínimo de error cuadrático segmentado; "
        f"mínimo R² postquiebre={MIN_POST_R2:.2f}; "
        f"mínimo de eventos por segmento={MIN_SEGMENT_EVENTS}; "
        "sin umbral mínimo total de eventos."
    )

    scale_definitions = (
        ("Regional (260 km)", "regional", "Catálogo regional", REGIONAL_RADIUS_KM),
        ("Local (50 km)", "local", "Catálogo local", LOCAL_RADIUS_KM),
    )
    completeness_tables = []
    annual_tables = []
    recurrence_tables = []
    corrected_fit_tables = []
    increment_tables = []
    raw_recurrence_tables = []
    raw_fit_tables = []

    for scale_name, scale_slug, catalog_name, radius in scale_definitions:
        scale_data = data[data["distancia_neiva_km"] <= radius].copy()
        if scale_data.empty:
            raise ValueError(
                f"No hay eventos en la escala {scale_name} (radio {radius:g} km)."
            )
        print(
            f"\nInicio del análisis de recurrencia: {scale_name}; "
            f"{len(scale_data):,} eventos."
        )
        completeness, fits = _make_completeness_tables(
            scale_data, years, sources, catalog_name, scale_name
        )
        completeness_tables.append(completeness)
        annual_path = tables_dir / (
            f"acumulado_anual_por_tramo_magnitud_{scale_slug}.csv"
        )
        _save_annual_cumulative_counts(
            scale_data,
            years,
            sources,
            annual_path,
            catalog_name,
            scale_name,
        )
        annual_tables.append(pd.read_csv(annual_path, encoding="utf-8-sig"))
        _save_completeness_figures(
            years,
            sources,
            fits,
            figures_dir,
            catalog_name,
            scale_slug,
        )
        recurrence, corrected_fits, increments = _make_recurrence_tables(
            scale_data,
            sources,
            completeness,
            first_year,
            last_complete_year,
            catalog_name,
            scale_name,
        )
        raw_recurrence, raw_fits = _make_uncorrected_recurrence(
            scale_data,
            sources,
            first_year,
            last_complete_year,
            catalog_name,
            scale_name,
        )
        recurrence_tables.append(recurrence)
        corrected_fit_tables.append(corrected_fits)
        increment_tables.append(increments)
        raw_recurrence_tables.append(raw_recurrence)
        raw_fit_tables.append(raw_fits)

    completeness = pd.concat(completeness_tables, ignore_index=True)
    annual = pd.concat(annual_tables, ignore_index=True)
    recurrence = pd.concat(recurrence_tables, ignore_index=True)
    gr_fits = pd.concat(corrected_fit_tables, ignore_index=True)
    increments = pd.concat(increment_tables, ignore_index=True)
    raw_recurrence = pd.concat(raw_recurrence_tables, ignore_index=True)
    raw_fits = pd.concat(raw_fit_tables, ignore_index=True)
    raw_fit_summary = raw_fits.rename(
        columns={
            "a": "a sin corrección",
            "b": "b sin corrección",
            "R2": "R2 sin corrección",
            "N intervalos ajuste": "N intervalos sin corrección",
        }
    )
    corrected_fit_summary = gr_fits.rename(
        columns={
            "a": "a corregido",
            "b": "b corregido",
            "R2": "R2 corregido",
            "N intervalos ajuste": "N intervalos corregidos",
        }
    )
    fit_comparison = raw_fit_summary.merge(
        corrected_fit_summary,
        on=["Escala", "Conjunto"],
        how="outer",
        validate="one_to_one",
    )

    completeness.to_csv(
        tables_dir / "completez_por_tramo_magnitud.csv",
        index=False,
        encoding="utf-8-sig",
    )
    annual.to_csv(
        tables_dir / "acumulado_anual_por_tramo_magnitud.csv",
        index=False,
        encoding="utf-8-sig",
    )
    recurrence.to_csv(
        tables_dir / "tabla_recurrencia_gr_por_tramo.csv",
        index=False,
        encoding="utf-8-sig",
    )
    gr_fits.to_csv(
        tables_dir / "ajuste_gutenberg_richter_sismofuentes.csv",
        index=False,
        encoding="utf-8-sig",
    )
    increments.to_csv(
        tables_dir / "conteos_incrementales_ventana_umbral.csv",
        index=False,
        encoding="utf-8-sig",
    )
    raw_recurrence.to_csv(
        tables_dir / "recurrencia_sin_correccion_por_tramo.csv",
        index=False,
        encoding="utf-8-sig",
    )
    fit_comparison.to_csv(
        tables_dir / "ajustes_crudos_y_corregidos_gutenberg_richter.csv",
        index=False,
        encoding="utf-8-sig",
    )
    _save_source_comparison_figures(
        recurrence, gr_fits, raw_recurrence, raw_fits, sources, figures_dir
    )
    _save_recurrence_plot(
        recurrence, gr_fits, raw_recurrence, figures_dir
    )

    print("\nCompletez por intervalos:")
    print(
        completeness[
            [
                "Escala",
                "Conjunto",
                "Intervalo Mw",
                "N total catálogo",
                "Año de inicio completo",
                "t completo (años)",
                "R2 postquiebre",
                "Estado",
            ]
        ].to_string(index=False)
    )
    print("\nAjustes Gutenberg--Richter:")
    print(gr_fits.to_string(index=False))
    print("\nAjustes sin corrección:")
    print(raw_fits.to_string(index=False))
    print(
        "\nEl análisis de recurrencia y la generación de todas sus figuras "
        "terminaron antes de ejecutar el análisis de atenuación."
    )
    print(f"\nTablas: {tables_dir}")
    print(f"Figuras: {figures_dir}")
    return {
        "completez": completeness,
        "recurrencia": recurrence,
        "ajustes": gr_fits,
        "recurrencia_sin_correccion": raw_recurrence,
        "ajustes_sin_correccion": raw_fits,
        "comparacion_ajustes": fit_comparison,
        "conteos_anuales": annual,
    }


if __name__ == "__main__":
    run_analysis()
