# -*- coding: utf-8 -*-
"""
Revisa formatos de fecha, separadores decimales y profundidades reportadas
en el catalogo consolidado (antes de la depuracion GK).
"""

from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DATOS = ROOT / "datos" / "datos depurados"

for radio in ["50KM", "260KM"]:
    archivo = DATOS / f"SISMO NEIVA-CONSOLIDADO-{radio}-HOMOGENEIZADO-MW.csv"
    print("=" * 70)
    print(f"RADIO: {radio}")
    print("=" * 70)

    # Leemos TODO como texto (dtype=str) para ver el formato crudo, sin que
    # pandas ya haya interpretado/limpiado los numeros por nosotros.
    crudo = pd.read_csv(archivo, encoding="utf-8-sig", dtype=str)

    print("\n--- Separadores decimales (columna 'magnitude', muestra) ---")
    for fuente in ["SGC", "USGS"]:
        subset = crudo[crudo["fuente"] == fuente]
        # Buscamos si hay comas en vez de puntos como separador decimal
        con_coma = subset["magnitude"].dropna().str.contains(",", na=False)
        print(f"{fuente}: {con_coma.sum()} valores con coma de {len(subset)} totales")
        print(f"  Ejemplos de magnitude: {subset['magnitude'].dropna().head(3).tolist()}")

    print("\n--- Separadores decimales (columna 'depth', muestra) ---")
    for fuente in ["SGC", "USGS"]:
        subset = crudo[crudo["fuente"] == fuente]
        con_coma = subset["depth"].dropna().str.contains(",", na=False)
        print(f"{fuente}: {con_coma.sum()} valores con coma de {len(subset)} totales")
        print(f"  Ejemplos de depth: {subset['depth'].dropna().head(3).tolist()}")

    # Ahora convertimos a numerico normal para revisar rangos de profundidad
    df = pd.read_csv(archivo, encoding="utf-8-sig")
    df["depth"] = pd.to_numeric(df["depth"], errors="coerce")

    print("\n--- Estadisticas de profundidad (depth) por fuente ---")
    for fuente in ["SGC", "USGS"]:
        subset = df[df["fuente"] == fuente]["depth"]
        print(f"{fuente}: N={len(subset)}  vacios={subset.isna().sum()}  "
              f"min={subset.min():.2f}  max={subset.max():.2f}  "
              f"negativos={(subset < 0).sum()}  ceros={(subset == 0).sum()}")

    print()