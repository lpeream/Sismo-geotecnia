# -*- coding: utf-8 -*-
"""
Busca posibles eventos duplicados entre SGC y USGS: el mismo sismo real
reportado por ambas fuentes con IDs distintos, detectado por cercania en
tiempo, ubicacion y magnitud.
"""

from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DATOS = ROOT / "datos" / "datos depurados"

# Umbrales de coincidencia
VENTANA_MINUTOS = 5
DISTANCIA_MAX_KM = 15
DIFERENCIA_MAG_MAX = 0.5


def distancia_km(lat1, lon1, lat2, lon2):
    radio_tierra = 6371.0
    lat1, lat2 = np.radians(lat1), np.radians(lat2)
    dlat = lat2 - lat1
    dlon = np.radians(lon2 - lon1)
    a = np.sin(dlat / 2) ** 2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2) ** 2
    return 2 * radio_tierra * np.arcsin(np.sqrt(a))


for radio in ["50KM", "260KM"]:
    archivo = DATOS / f"SISMO NEIVA-CONSOLIDADO-{radio}-HOMOGENEIZADO-MW.csv"
    df = pd.read_csv(archivo, encoding="utf-8-sig")
    df["fecha_utc"] = pd.to_datetime(df["fecha_utc"], errors="coerce", format="ISO8601")
    df["magnitude_mw"] = pd.to_numeric(df["magnitude_mw"], errors="coerce")
    df["latitude"] = pd.to_numeric(df["latitude"], errors="coerce")
    df["longitude"] = pd.to_numeric(df["longitude"], errors="coerce")

    sgc = df[df["fuente"] == "SGC"].dropna(subset=["fecha_utc", "latitude", "longitude", "magnitude_mw"])
    usgs = df[df["fuente"] == "USGS"].dropna(subset=["fecha_utc", "latitude", "longitude", "magnitude_mw"])

    print("=" * 60)
    print(f"RADIO: {radio}  (SGC: {len(sgc)}, USGS: {len(usgs)})")
    print("=" * 60)

    posibles_duplicados = []
    for _, evento_usgs in usgs.iterrows():
        ventana = sgc[
            (sgc["fecha_utc"] >= evento_usgs["fecha_utc"] - pd.Timedelta(minutes=VENTANA_MINUTOS))
            & (sgc["fecha_utc"] <= evento_usgs["fecha_utc"] + pd.Timedelta(minutes=VENTANA_MINUTOS))
        ]
        if ventana.empty:
            continue
        distancias = distancia_km(
            evento_usgs["latitude"], evento_usgs["longitude"],
            ventana["latitude"].to_numpy(), ventana["longitude"].to_numpy(),
        )
        cercanos = ventana[distancias <= DISTANCIA_MAX_KM].copy()
        if cercanos.empty:
            continue
        cercanos["dif_magnitud"] = (cercanos["magnitude_mw"] - evento_usgs["magnitude_mw"]).abs()
        cercanos = cercanos[cercanos["dif_magnitud"] <= DIFERENCIA_MAG_MAX]
        for _, evento_sgc in cercanos.iterrows():
            posibles_duplicados.append({
                "id_usgs": evento_usgs["id"], "fecha_usgs": evento_usgs["fecha_utc"],
                "mag_usgs": evento_usgs["magnitude_mw"],
                "id_sgc": evento_sgc["id"], "fecha_sgc": evento_sgc["fecha_utc"],
                "mag_sgc": evento_sgc["magnitude_mw"],
            })

    print(f"Posibles duplicados encontrados: {len(posibles_duplicados)}")
    if posibles_duplicados:
        tabla = pd.DataFrame(posibles_duplicados)
        print(tabla.to_string())
        tabla.to_csv(DATOS / f"posibles_duplicados_{radio.lower()}.csv", index=False, encoding="utf-8-sig")
        print(f"-> Guardado en posibles_duplicados_{radio.lower()}.csv")
    print()