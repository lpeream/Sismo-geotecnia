from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


DELIVERY_ROOT = Path(__file__).resolve().parents[1]
HAZARD_CSV = (
    DELIVERY_ROOT
    / "datos"
    / "datos estadistica"
    / "integracion paso4"
    / "curva_amenaza_neiva_roca_paso4.csv"
)
OUTPUT_PNG = (
    DELIVERY_ROOT
    / "figuras"
    / "amenaza_sismica_neiva_figura8_local_aritmetica.png"
)

PGA_MAX = 20.0
RATE_MAX = 20.0
SOURCES = (
    {
        "label": "Sismofuente Cortical",
        "column": "lambda_Cortical_por_anio",
        "color": "#C00000",
        "marker": "D",
    },
    {
        "label": "Subducción Intraplaca",
        "column": "lambda_Subducción intraplaca_por_anio",
        "color": "#002060",
        "marker": "^",
    },
    {
        "label": "Subducción Interplaca",
        "column": "lambda_Subducción interplaca_por_anio",
        "color": "#7030A0",
        "marker": "s",
    },
    {
        "label": "Amenaza Total en Roca",
        "column": "lambda_total_por_anio",
        "color": "black",
        "marker": None,
    },
)


def load_hazard_data() -> tuple[np.ndarray, dict[str, np.ndarray]]:
    if not HAZARD_CSV.is_file():
        raise FileNotFoundError(
            f"No se encontró la curva calculada de amenaza: {HAZARD_CSV}"
        )

    required_columns = {
        "PGA_roca_cm_s2",
        *(source["column"] for source in SOURCES),
    }
    with HAZARD_CSV.open(encoding="utf-8-sig", newline="") as file:
        reader = csv.DictReader(file)
        if reader.fieldnames is None:
            raise ValueError(f"El archivo no contiene encabezados: {HAZARD_CSV}")
        missing_columns = required_columns.difference(reader.fieldnames)
        if missing_columns:
            raise ValueError(
                "Faltan columnas en la curva de amenaza: "
                + ", ".join(sorted(missing_columns))
            )
        rows = list(reader)

    if not rows:
        raise ValueError(f"La curva de amenaza está vacía: {HAZARD_CSV}")

    try:
        pga = np.array(
            [float(row["PGA_roca_cm_s2"]) for row in rows],
            dtype=float,
        )
        rates = {
            source["column"]: np.array(
                [float(row[source["column"]]) for row in rows],
                dtype=float,
            )
            for source in SOURCES
        }
    except (TypeError, ValueError) as error:
        raise ValueError(
            f"La curva de amenaza contiene valores no numéricos: {HAZARD_CSV}"
        ) from error

    visible = pga <= PGA_MAX
    pga = pga[visible]
    rates = {column: values[visible] for column, values in rates.items()}
    if pga.size == 0 or np.any(np.diff(pga) <= 0.0):
        raise ValueError("Los valores de PGA deben ser estrictamente crecientes.")
    if any(not np.isfinite(values).all() for values in rates.values()):
        raise ValueError("La curva de amenaza contiene tasas no finitas.")
    if any(np.any(values < 0.0) for values in rates.values()):
        raise ValueError("La curva de amenaza contiene tasas negativas.")

    source_total = sum(
        rates[source["column"]]
        for source in SOURCES
        if source["column"] != "lambda_total_por_anio"
    )
    if not np.allclose(source_total, rates["lambda_total_por_anio"]):
        raise ArithmeticError(
            "La tasa total no coincide con la suma de las tres sismofuentes."
        )
    return pga, rates


def main() -> None:
    pga, rates = load_hazard_data()
    figure, axis = plt.subplots(figsize=(10, 6.2), dpi=300)
    figure.patch.set_facecolor("white")
    axis.set_facecolor("white")

    for source in SOURCES:
        rate = rates[source["column"]]
        axis.plot(
            pga,
            rate,
            color=source["color"],
            linestyle="-" if source["marker"] is None else "--",
            linewidth=2.5 if source["marker"] is None else 1.6,
            label=source["label"],
        )
        if source["marker"] is not None:
            axis.plot(
                pga,
                rate,
                linestyle="none",
                marker=source["marker"],
                markevery=16,
                markersize=4,
                markerfacecolor="none",
                markeredgecolor=source["color"],
                markeredgewidth=1.1,
            )

    total_rate_at_zero = float(rates["lambda_total_por_anio"][0])
    axis.text(
        0.98,
        0.97,
        rf"$\lambda_{{\mathrm{{total}}}}(0)={total_rate_at_zero:.2f}"
        r"\ \mathrm{año}^{-1}$ (fuera de escala)",
        transform=axis.transAxes,
        horizontalalignment="right",
        verticalalignment="top",
        color="black",
        fontsize=8,
    )
    axis.set_xlim(0.0, PGA_MAX)
    axis.set_ylim(0.0, RATE_MAX)
    axis.set_xlabel("PGA [cm/s²]", fontsize=11, fontweight="bold")
    axis.set_ylabel(
        r"$\lambda$",
        fontsize=13,
        fontweight="bold",
        rotation=0,
        labelpad=14,
    )
    axis.set_title(
        "Amenaza sísmica integrada en roca - Neiva",
        fontsize=14,
        fontweight="bold",
    )
    axis.grid(color="#D9D9D9", linestyle=":", linewidth=0.8)
    axis.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, -0.12),
        ncol=2,
        frameon=False,
    )
    figure.tight_layout()
    figure.savefig(OUTPUT_PNG, dpi=300, bbox_inches="tight")
    plt.close(figure)
    print(f"Figura guardada en: {OUTPUT_PNG}")


if __name__ == "__main__":
    main()
