# -*- coding: utf-8 -*-
"""
Consolida los catalogos de SGC y USGS y homogeniza las magnitudes a Mw.

Las magnitudes originales se conservan. Para cada evento se agrega un Mw
directo o estimado, la correlacion aplicada, su rango y una bandera de calidad.
Los tipos de magnitud ambiguos no se convierten automaticamente.
"""

import pandas as pd

ORIGINALES = "../datos/datos originales"
SALIDA = "../datos/datos depurados"

COLUMNAS_COMUNES = [
    "fuente", "id", "fecha_utc", "fecha_local", "latitude", "longitude",
    "depth", "magnitude", "mag_type", "place", "status", "rms", "gap", "nst",
]


def convertir_a_mw(magnitude, mag_type):
    """Devuelve Mw, relacion, sigma y estado de la conversion.

    Las relaciones corresponden a las ecuaciones incluidas en main.tex:
    Scordilis/Bormann-Saul para mb, Ms y MJMA, y Senel et al. para ML y Md.
    """
    try:
        valor = float(magnitude)
    except (TypeError, ValueError):
        return (pd.NA, "sin magnitud", pd.NA, "sin_magnitud")

    tipo = "" if pd.isna(mag_type) else str(mag_type).strip().lower()
    tipo_base = tipo.replace("_", "").replace("-", "")

    # Estas escalas ya son estimaciones de magnitud de momento.
    if tipo in {"mw", "mww", "mwc", "mwb", "mwr", "mwp"}:
        return (valor, "directa: " + str(mag_type), 0.0, "directa")
    if tipo in {"mw(mb)", "mwb(mb)"}:
        return (valor, "directa reportada: " + str(mag_type), 0.0, "directa")

    # Scordilis (2005), reproducido por Bormann y Saul (2005).
    if tipo in {"mb", "mbasic"} or tipo_base == "mb":
        mw = 0.85 * valor + 1.02
        estado = "dentro_de_rango" if 3.5 <= valor <= 6.2 else "extrapolada"
        return (mw, "Scordilis: Mw=0.85 mb+1.02 (3.5<=mb<=6.2)", 0.29, estado)

    if tipo in {"ms", "msk"}:
        if valor <= 6.1:
            mw = 0.65 * valor + 2.20
            sigma = 0.17
            rango = "3.0<=Ms<=6.1"
        else:
            mw = 1.00 * valor - 0.02
            sigma = 0.21
            rango = "6.2<=Ms<=8.0"
        estado = "dentro_de_rango" if 3.0 <= valor <= 8.0 else "extrapolada"
        return (mw, f"Scordilis: {rango}", sigma, estado)

    # Senel et al. (2016): unica relacion explicita del corpus para ML/Md.
    es_local = (
        tipo in {"ml", "mlr", "mlv"}
        or tipo.startswith("mlr")
        or "mlr" in tipo
    )
    if es_local:
        mw = 0.8095 * valor + 1.3003
        estado = "dentro_de_rango" if 3.3 <= valor <= 6.6 else "extrapolada"
        return (mw, "Senel et al.: Mw=0.8095 ML+1.3003 (3.3<=ML<=6.6)", pd.NA, estado)

    if tipo in {"md", "mc"}:
        mw = 0.7947 * valor + 1.3420
        estado = "dentro_de_rango" if 3.5 <= valor <= 7.4 else "extrapolada"
        return (mw, "Senel et al.: Mw=0.7947 Md+1.3420 (3.5<=Md<=7.4)", pd.NA, estado)

    return (pd.NA, "sin correlacion verificable", pd.NA, "no_convertida")


def agregar_mw(df):
    conversiones = [
        convertir_a_mw(magnitude, mag_type)
        for magnitude, mag_type in zip(df["magnitude"], df["mag_type"])
    ]
    df[["magnitude_mw", "conversion_relation", "conversion_sigma",
        "conversion_status"]] = pd.DataFrame(conversiones, index=df.index)
    return df


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


for radio in ["50KM", "260KM"]:
    sgc = cargar_sgc(radio)
    usgs = cargar_usgs(radio)
    consolidado = pd.concat([sgc, usgs], ignore_index=True, sort=False)
    consolidado = consolidado.sort_values("fecha_utc").reset_index(drop=True)
    consolidado = agregar_mw(consolidado)

    archivo_salida = f"{SALIDA}/SISMO NEIVA-CONSOLIDADO-{radio}-HOMOGENEIZADO-MW.csv"
    consolidado.to_csv(archivo_salida, index=False, encoding="utf-8-sig")

    print(f"{radio}: SGC={len(sgc)}  USGS={len(usgs)}  consolidado={len(consolidado)}")
    print(f"  -> {archivo_salida}")
    print("  estados:")
    print(consolidado["conversion_status"].value_counts(dropna=False).to_string())
