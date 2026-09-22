# -*- coding: utf-8 -*-
"""
Genera el mapa de la zona de estudio (Neiva, Huila) con el punto de
referencia, los radios de busqueda de 50 y 260 km, y los epicentros
del catalogo consolidado (SGC + USGS, 260 km) coloreados por magnitud
y diferenciados por fuente.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

REF_LAT = 2.9262767
REF_LON = -75.2892111
KM_POR_GRADO_LAT = 111.32

DATOS = "../datos/datos depurados/SISMO NEIVA-CONSOLIDADO-260KM-SIN-FILTRAR.csv"
SALIDA = "../figuras/mapa_zona_estudio.png"


def circulo(lat0, lon0, radio_km, n=200):
    km_por_grado_lon = KM_POR_GRADO_LAT * np.cos(np.radians(lat0))
    t = np.linspace(0, 2 * np.pi, n)
    lat = lat0 + (radio_km / KM_POR_GRADO_LAT) * np.sin(t)
    lon = lon0 + (radio_km / km_por_grado_lon) * np.cos(t)
    return lon, lat


df = pd.read_csv(DATOS)
df["magnitude"] = pd.to_numeric(df["magnitude"], errors="coerce")
df = df.dropna(subset=["magnitude", "latitude", "longitude"])

sgc = df[df["fuente"] == "SGC"].sort_values("magnitude")
usgs = df[df["fuente"] == "USGS"].sort_values("magnitude")

fig, ax = plt.subplots(figsize=(8.5, 8))

vmin, vmax = df["magnitude"].min(), df["magnitude"].max()

sc = ax.scatter(
    sgc["longitude"], sgc["latitude"],
    c=sgc["magnitude"], s=8 + sgc["magnitude"] ** 2.2,
    cmap="inferno_r", alpha=0.45, linewidths=0,
    vmin=vmin, vmax=vmax, zorder=2, marker="o",
    label=f"SGC (n={len(sgc)})",
)

ax.scatter(
    usgs["longitude"], usgs["latitude"],
    c=usgs["magnitude"], s=25 + usgs["magnitude"] ** 2.6,
    cmap="inferno_r", alpha=0.95, linewidths=0.6, edgecolors="black",
    vmin=vmin, vmax=vmax, zorder=3, marker="^",
    label=f"USGS (n={len(usgs)})",
)

for radio, estilo in [(50, "-"), (260, "--")]:
    lon_c, lat_c = circulo(REF_LAT, REF_LON, radio)
    ax.plot(lon_c, lat_c, estilo, color="black", linewidth=1.2, zorder=4)
    ax.annotate(
        f"{radio} km",
        (lon_c[len(lon_c) // 8], lat_c[len(lat_c) // 8]),
        fontsize=9, color="black", zorder=5,
    )

ax.scatter(
    [REF_LON], [REF_LAT], marker="*", s=280, color="#1f77ff",
    edgecolor="black", linewidth=0.8, zorder=6,
    label="Punto de referencia\n(Parque Central Santander, Neiva)",
)

ax.set_xlabel("Longitud (°)")
ax.set_ylabel("Latitud (°)")
ax.set_title(
    "Zona de estudio: Neiva, Huila\n"
    "Catálogo consolidado SGC + USGS (260 km), por magnitud y fuente",
    fontsize=12,
)
ax.set_aspect(1 / np.cos(np.radians(REF_LAT)))
ax.grid(alpha=0.25, zorder=0)
ax.legend(loc="upper left", fontsize=8, framealpha=0.9, markerscale=0.9)

cbar = fig.colorbar(sc, ax=ax, shrink=0.8, pad=0.02)
cbar.set_label("Magnitud (escala original de cada fuente)")

fig.tight_layout()
fig.savefig(SALIDA, dpi=200)
print(f"Guardado: {SALIDA}")
print(f"Eventos SGC: {len(sgc)}  |  Eventos USGS: {len(usgs)}")
