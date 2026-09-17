# -*- coding: utf-8 -*-
"""
Consolida (concatena) los catalogos de SGC y USGS para un mismo radio,
en una sola tabla con columnas comunes. NO homogeniza magnitudes a Mw,
NO filtra por magnitud minima y NO elimina replicas: solo unifica las
dos fuentes conservando toda la informacion original de cada una.
"""

import pandas as pd

ORIGINALES = "../datos/datos originales"
SALIDA = "../datos/datos depurados"

COLUMNAS_COMUNES = [
    "fuente", "id", "fecha_utc", "fecha_local", "latitude", "longitude",
    "depth", "magnitude", "mag_type", "place", "status", "rms", "gap", "nst",
]


def cargar_sgc(radio):
    df = pd.read_csv(f"{ORIGINALES}/SISMO NEIVA-SGC-{radio}-FULL.csv")
    out = pd.DataFrame({
        "fuente": "SGC",
        "id": df["id"],
        "fecha_utc": df["utc_time"],
        "fecha_local": df["local_time"],
        "latitude": df["latitude"],
        "longitude": df["longitude"],
        "depth": df["depth"],
        "magnitude": df["magnitude"],
        "mag_type": df["mag_type"],
        "place": df["place"],
        "status": df["status"],
        "rms": df["rms"],
        "gap": df["gap"],
        "nst": df["nst"],
    })
    # columnas propias del SGC, se conservan sin perder informacion
    out["agency"] = df["agency"]
    out["magnitude_error"] = df["magnitude_error"]
    out["latitude_error"] = df["latitude_error"]
    out["longitude_error"] = df["longitude_error"]
    out["depth_error"] = df["depth_error"]
    return out


def cargar_usgs(radio):
    df = pd.read_csv(f"{ORIGINALES}/SISMO NEIVA-USGS-{radio}.csv")
    out = pd.DataFrame({
        "fuente": "USGS",
        "id": df["id"],
        "fecha_utc": df["time"],
        "fecha_local": pd.NA,
        "latitude": df["latitude"],
        "longitude": df["longitude"],
        "depth": df["depth"],
        "magnitude": df["mag"],
        "mag_type": df["magType"],
        "place": df["place"],
        "status": df["status"],
        "rms": df["rms"],
        "gap": df["gap"],
        "nst": df["nst"],
    })
    # columnas propias del USGS, se conservan sin perder informacion
    out["net"] = df["net"]
    out["type"] = df["type"]
    out["horizontalError"] = df["horizontalError"]
    out["depthError"] = df["depthError"]
    out["magError"] = df["magError"]
    out["magNst"] = df["magNst"]
    out["locationSource"] = df["locationSource"]
    out["magSource"] = df["magSource"]
    out["dmin"] = df["dmin"]
    return out


for radio in ["50KM", "200KM"]:
    sgc = cargar_sgc(radio)
    usgs = cargar_usgs(radio)
    consolidado = pd.concat([sgc, usgs], ignore_index=True, sort=False)
    consolidado = consolidado.sort_values("fecha_utc").reset_index(drop=True)

    archivo_salida = f"{SALIDA}/SISMO NEIVA-CONSOLIDADO-{radio}-SIN-FILTRAR.csv"
    consolidado.to_csv(archivo_salida, index=False, encoding="utf-8-sig")

    print(f"{radio}: SGC={len(sgc)}  USGS={len(usgs)}  consolidado={len(consolidado)}")
    print(f"  -> {archivo_salida}")
