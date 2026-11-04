from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

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
    [3.0, 3.3, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5]
)
MIN_SEGMENT_YEARS = 10
MIN_SEGMENT_EVENTS = 3
MIN_EVENTS_FOR_BREAK = 10
MIN_POST_R2 = 0.90
MIN_FIT_BINS = 3


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

    if n_events < MIN_EVENTS_FOR_BREAK:
        return unestimated("No estimable: menos de 10 eventos")

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
        return unestimated("No estimable: no cumple mínimos por segmento")

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
) -> tuple[pd.DataFrame, dict[tuple[str, str], CompletenessFit]]:
    rows = []
    fits: dict[tuple[str, str], CompletenessFit] = {}
    groups = [("Catálogo regional", data)] + [
        (source, data[data["source"] == source]) for source in group_names
    ]

    for group_name, group in groups:
        for lower, upper, label in _magnitude_bins():
            band = group[group["Mw"].between(lower, upper, inclusive="left")]
            fit = _fit_completeness(band, years)
            fits[(group_name, label)] = fit
            start_year = fit.start_year
            rows.append(
                {
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
    data: pd.DataFrame, years: np.ndarray, groups: list[str], output_path: Path
) -> None:
    rows = []
    all_groups = [("Catálogo regional", data)] + [
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
) -> None:
    all_groups = ["Catálogo regional", *groups]
    for group_name in all_groups:
        fig, axes = plt.subplots(5, 2, figsize=(12, 16), sharex=True)
        for ax, (_, _, label) in zip(axes.flat, _magnitude_bins()):
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
        axes[-1, 0].set_xlabel("Año")
        axes[-1, 1].set_xlabel("Año")
        axes[0, 0].legend(fontsize=7, loc="upper left")
        fig.suptitle(
            f"Completez por intervalos de magnitud — {group_name}",
            y=0.995,
        )
        fig.tight_layout(rect=(0, 0, 1, 0.99))
        slug = (
            "catalogo_regional"
            if group_name == "Catálogo regional"
            else group_name.lower()
            .replace(" ", "_")
            .replace("ó", "o")
        )
        fig.savefig(
            figures_dir / f"completez_acumulada_{slug}.png",
            dpi=220,
        )
        plt.close(fig)


def _make_recurrence_tables(
    data: pd.DataFrame,
    group_names: list[str],
    completeness: pd.DataFrame,
    last_complete_year: int,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    recurrence_rows = []
    fit_rows = []
    increment_rows = []
    fit_data: dict[str, list[dict[str, float]]] = {}
    groups = [("Catálogo regional", data)] + [
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
            }
            if pd.notna(start_year_value):
                start_year = int(start_year_value)
                duration = last_complete_year - start_year + 1
                observed = group[group["year"].between(start_year, last_complete_year)]
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
                        "log10(lambda)": (
                            float(np.log10(rate)) if rate > 0 else np.nan
                        ),
                    }
                )
                if rate > 0:
                    valid_points.append(
                        {"Mw mínimo": lower, "log10(lambda)": float(np.log10(rate))}
                    )
            recurrence_rows.append(row)

        if len(valid_points) >= MIN_FIT_BINS:
            x = np.array([point["Mw mínimo"] for point in valid_points])
            y = np.array([point["log10(lambda)"] for point in valid_points])
            slope, intercept = np.polyfit(x, y, 1)
            predicted = intercept + slope * x
            r_squared = _r_squared(y, predicted)
            fit = {
                "Conjunto": group_name,
                "a": float(intercept),
                "b": float(-slope),
                "R2": r_squared,
                "N intervalos ajuste": len(valid_points),
            }
            fit_rows.append(fit)
            fit_data[group_name] = valid_points
        else:
            fit = {
                "Conjunto": group_name,
                "a": np.nan,
                "b": np.nan,
                "R2": np.nan,
                "N intervalos ajuste": len(valid_points),
            }
            fit_rows.append(fit)
            fit_data[group_name] = []

    recurrence = pd.DataFrame(recurrence_rows)
    fits = pd.DataFrame(fit_rows)
    increments = pd.DataFrame(increment_rows)
    for fit in fit_rows:
        selected = recurrence["Conjunto"] == fit["Conjunto"]
        recurrence.loc[selected, "a"] = fit["a"]
        recurrence.loc[selected, "b"] = fit["b"]
        recurrence.loc[selected, "R2 ajuste"] = fit["R2"]
    return recurrence, fits, increments


def _save_recurrence_plot(
    recurrence: pd.DataFrame,
    fits: pd.DataFrame,
    figures_dir: Path,
    last_complete_year: int,
) -> None:
    fig, ax = plt.subplots(figsize=(8, 5.5))
    for fit in fits.itertuples(index=False):
        subset = recurrence[
            (recurrence["Conjunto"] == fit.Conjunto)
            & recurrence["log10(lambda)"].notna()
        ]
        if subset.empty:
            continue
        ax.plot(
            subset["Mw mínimo"],
            subset["lambda acumulada (año^-1)"],
            "o",
            label=f"{fit.Conjunto} (b={fit.b:.2f})",
        )
        if pd.notna(fit.a) and pd.notna(fit.b):
            x_fit = np.linspace(
                float(subset["Mw mínimo"].min()),
                float(subset["Mw mínimo"].max()),
                100,
            )
            ax.plot(x_fit, 10 ** (fit.a - fit.b * x_fit), "--")
    ax.set_yscale("log")
    ax.set_xlabel(r"Umbral de magnitud $M_w$")
    ax.set_ylabel(
        r"Tasa anual acumulada $\lambda(M_w\geq M)$ (eventos/año)"
    )
    ax.set_title(
        f"Recurrencia Gutenberg--Richter con completez por tramos "
        f"(hasta {last_complete_year})"
    )
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(
        figures_dir / "curvas_recurrencia_sismofuentes.png", dpi=250
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
    required = {"Mw", "hora", "fuente"}
    missing_columns = required.difference(catalog.columns)
    if missing_columns:
        raise ValueError(
            "Faltan columnas requeridas: " + ", ".join(sorted(missing_columns))
        )
    catalog["Mw"] = pd.to_numeric(catalog["Mw"], errors="coerce")
    catalog["date_utc"] = pd.to_datetime(
        catalog["hora"], format="mixed", errors="coerce", utc=True
    )
    if catalog[["Mw", "date_utc", "fuente"]].isna().any().any():
        raise ValueError("El catálogo contiene magnitudes, fechas o fuentes vacías.")
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
        f"mínimo de eventos por intervalo={MIN_EVENTS_FOR_BREAK}."
    )

    completeness, fits = _make_completeness_tables(data, years, sources)
    completeness.to_csv(
        tables_dir / "completez_por_tramo_magnitud.csv",
        index=False,
        encoding="utf-8-sig",
    )
    _save_annual_cumulative_counts(
        data,
        years,
        sources,
        tables_dir / "acumulado_anual_por_tramo_magnitud.csv",
    )
    _save_completeness_figures(years, sources, fits, figures_dir)

    recurrence, gr_fits, increments = _make_recurrence_tables(
        data, sources, completeness, last_complete_year
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
    _save_recurrence_plot(
        recurrence, gr_fits, figures_dir, last_complete_year
    )

    print("\nCompletez por intervalos:")
    print(
        completeness[
            [
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
    print(f"\nTablas: {tables_dir}")
    print(f"Figuras: {figures_dir}")
    return {
        "completez": completeness,
        "recurrencia": recurrence,
        "ajustes": gr_fits,
    }


if __name__ == "__main__":
    run_analysis()
