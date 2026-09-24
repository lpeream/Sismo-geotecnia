# -*- coding: utf-8 -*-
"""
Revisa manualmente los registros marcados como no_convertida o extrapolada
en el catalogo homogeneizado, y genera un resumen para inspeccion.
"""

import pandas as pd

DATOS = "../datos/datos depurados"

pd.set_option("display.max_columns", None)
pd.set_option("display.width", 200)

for radio in ["50KM", "260KM"]:
    archivo = f"{DATOS}/SISMO NEIVA-CONSOLIDADO-{radio}-HOMOGENEIZADO-MW.csv"
    df = pd.read_csv(archivo)

    print(f"\n{'='*60}")
    print(f"RADIO: {radio}")
    print(f"{'='*60}")
    print(f"Total de registros: {len(df)}")
    print("\nDistribución de conversion_status:")
    print(df["conversion_status"].value_counts(dropna=False))

    no_convertidas = df[df["conversion_status"] == "no_convertida"]
    extrapoladas = df[df["conversion_status"] == "extrapolada"]

    print(f"\nRegistros no_convertida: {len(no_convertidas)}")
    if len(no_convertidas) > 0:
        print("Tipos de magnitud (mag_type) involucrados:")
        print(no_convertidas["mag_type"].value_counts(dropna=False))

    print(f"\nRegistros extrapolada: {len(extrapoladas)}")
    if len(extrapoladas) > 0:
        print("Tipos de magnitud (mag_type) involucrados:")
        print(extrapoladas["mag_type"].value_counts(dropna=False))
        print("\nRango de magnitudes originales extrapoladas:")
        print(extrapoladas.groupby("mag_type")["magnitude"].agg(["min", "max", "count"]))

    columnas_revision = [
        "fuente", "id", "fecha_utc", "magnitude", "mag_type",
        "magnitude_mw", "conversion_relation", "conversion_status", "place",
    ]
    salida = f"{DATOS}/SISMO NEIVA-REVISION-NOCONV-EXTRAP-{radio}.csv"
    revision = pd.concat([no_convertidas, extrapoladas])[columnas_revision]
    revision.to_csv(salida, index=False, encoding="utf-8-sig")
    print(f"\n  -> Exportado para revisión: {salida}")