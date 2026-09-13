# -*- coding: utf-8 -*-
"""
Genera el mapa de la zona de estudio (Neiva, Huila) con el punto de
referencia, los radios de busqueda de 50 y 200 km, y los epicentros
del catalogo SGC (200 km) coloreados por magnitud.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

REF_LAT = 2.9262767
REF_LON = -75.2892111
KM_POR_GRADO_LAT = 111.32

DATOS = "../datos/datos originales/SISMO NEIVA-SGC-200KM-FULL.csv"
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

fig, ax = plt.subplots(figsize=(7, 7))

orden = df["magnitude"].argsort()
sc = ax.scatter(
    df["longitude"].values[orden],
    df["latitude"].values[orden],
    c=df["magnitude"].values[orden],
    s=8 + df["magnitude"].values[orden] ** 2.2,
    cmap="inferno_r",
    alpha=0.55,
    linewidths=0,
    vmin=df["magnitude"].min(),
    vmax=df["magnitude"].max(),
    zorder=2,
)

for radio, estilo in [(50, "-"), (200, "--")]:
    lon_c, lat_c = circulo(REF_LAT, REF_LON, radio)
    ax.plot(lon_c, lat_c, estilo, color="black", linewidth=1.2, zorder=3)
    ax.annotate(
        f"{radio} km",
        (lon_c[len(lon_c) // 8], lat_c[len(lat_c) // 8]),
        fontsize=9,
        color="black",
        zorder=4,
    )

ax.scatter(
    [REF_LON], [REF_LAT], marker="*", s=260, color="#1f77ff",
    edgecolor="black", linewidth=0.8, zorder=5,
    label="Punto de referencia\n(Parque Central Santander, Neiva)",
)

ax.set_xlabel("Longitud (°)")
ax.set_ylabel("Latitud (°)")
ax.set_title("Zona de estudio: Neiva, Huila\nCatálogo SGC dentro de 200 km, coloreado por magnitud")
ax.set_aspect(1 / np.cos(np.radians(REF_LAT)))
ax.grid(alpha=0.25, zorder=0)
ax.legend(loc="upper left", fontsize=8, framealpha=0.9)

cbar = fig.colorbar(sc, ax=ax, shrink=0.8, pad=0.02)
cbar.set_label("Magnitud")

fig.tight_layout()
fig.savefig(SALIDA, dpi=200)
print(f"Guardado: {SALIDA}")
print(f"Eventos graficados: {len(df)}")
