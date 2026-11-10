from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


DELIVERY_ROOT = Path(__file__).resolve().parents[1]
OUTPUT_PNG = (
    DELIVERY_ROOT
    / "figuras"
    / "amenaza_exponencial_neiva_regional.png"
)

PGA_MAX_CM_S2 = 40.0
PGA_GRID = np.linspace(0.0, PGA_MAX_CM_S2, 401)
MARKER_PGA = np.arange(2.0, PGA_MAX_CM_S2 + 0.1, 4.0)

MODELS = (
    {
        "label": "Sismofuente Cortical",
        "amplitude": 24.850,
        "decay": 0.085,
        "r_squared": 0.9854,
        "color": "#C00000",
        "marker": "D",
    },
    {
        "label": "Subducción Intraplaca",
        "amplitude": 16.400,
        "decay": 0.112,
        "r_squared": 0.9912,
        "color": "#002060",
        "marker": "^",
    },
    {
        "label": "Subducción Interplaca",
        "amplitude": 6.250,
        "decay": 0.245,
        "r_squared": 0.9788,
        "color": "#7030A0",
        "marker": "s",
    },
    {
        "label": "Amenaza Total en Roca",
        "amplitude": 46.276,
        "decay": 0.101,
        "r_squared": 0.9815,
        "color": "black",
        "marker": None,
    },
)


def main() -> None:
    figure, axis = plt.subplots(figsize=(10, 6.2), dpi=300)
    figure.patch.set_facecolor("white")
    axis.set_facecolor("white")

    for model in MODELS:
        fitted_rate = model["amplitude"] * np.exp(
            -model["decay"] * PGA_GRID
        )
        axis.plot(
            PGA_GRID,
            fitted_rate,
            color=model["color"],
            linestyle="-" if model["marker"] is None else "--",
            linewidth=2.5 if model["marker"] is None else 1.6,
            label=model["label"],
        )
        if model["marker"] is not None:
            marker_rate = model["amplitude"] * np.exp(
                -model["decay"] * MARKER_PGA
            )
            axis.plot(
                MARKER_PGA,
                marker_rate,
                linestyle="none",
                marker=model["marker"],
                markersize=5,
                markerfacecolor="none",
                markeredgecolor=model["color"],
                markeredgewidth=1.1,
            )

    for index, model in enumerate(MODELS):
        annotation = (
            f"{model['label']}: "
            rf"$y={model['amplitude']:.3f}\,e^{{-{model['decay']:.3f}x}}$, "
            rf"$R^2={model['r_squared']:.4f}$"
        )
        axis.text(
            0.98,
            0.97 - index * 0.075,
            annotation,
            transform=axis.transAxes,
            horizontalalignment="right",
            verticalalignment="top",
            color=model["color"],
            fontsize=8,
        )

    axis.set_xlim(0.0, PGA_MAX_CM_S2)
    axis.set_ylim(0.0, 40.0)
    axis.set_xlabel("PGA [cm/s²]", fontsize=11, fontweight="bold")
    axis.set_ylabel(
        r"$\lambda$",
        fontsize=13,
        fontweight="bold",
        rotation=0,
        labelpad=14,
    )
    axis.set_title(
        "Curvas de amenaza sísmica probabilística en roca base - Neiva",
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
