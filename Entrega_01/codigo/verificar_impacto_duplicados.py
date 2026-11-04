# -*- coding: utf-8 -*-
"""
Verifica el impacto de eliminar los posibles duplicados SGC/USGS sobre el
ajuste de Gutenberg-Richter en el catalogo de 260 km.
"""

from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DATOS = ROOT / "datos" / "datos depurados"
ANCHO_BIN = 0.1


def ajustar_gr(magnitudes):
    minimo = np.floor(magnitudes.min() / ANCHO_BIN) * ANCHO_BIN
    maximo = np.ceil(magnitudes.max() / ANCHO_BIN) * ANCHO_BIN + ANCHO_BIN
    bordes = np.arange(minimo, maximo + ANCHO_BIN / 2, ANCHO_BIN)
    frecuencias, _ = np.histogram(magnitudes, bins=bordes)
    centros = bordes[:-1] + ANCHO_BIN / 2
    mc = float(centros[int(np.argmax(frecuencias))])
    valores = np.sort(magnitudes[magnitudes >= mc])
    magnitudes_unicas = np.unique(valores)
    acumuladas = np.array([(valores >= v).sum() for v in magnitudes_unicas])
    mascara = acumuladas >= 2
    x = magnitudes_unicas[mascara]
    y = np.log10(acumuladas[mascara])
    pendiente, intercepto = np.polyfit(x, y, 1)
    prediccion = intercepto + pendiente * x
    r2 = 1 - np.sum((y - prediccion) ** 2) / np.sum((y - y.mean()) ** 2)
    return {"mc": mc, "b": -pendiente, "r2": r2, "n_total": len(magnitudes), "n_ajuste": len(valores)}


duplicados = pd.read_csv(DATOS / "posibles_duplicados_260km.csv", encoding="utf-8-sig")
ids_usgs_duplicados = set(duplicados["id_usgs"].unique())
print(f"IDs unicos de USGS marcados como duplicados: {len(ids_usgs_duplicados)}")

archivo = DATOS / "SISMO NEIVA-CONSOLIDADO-260KM-HOMOGENEIZADO-MW-DEPURADO-GK.csv"
df = pd.read_csv(archivo, encoding="utf-8-sig")
df["magnitude_mw"] = pd.to_numeric(df["magnitude_mw"], errors="coerce")

print(f"\nCatalogo original (post Gardner-Knopoff): N={len(df)}")

antes = ajustar_gr(df["magnitude_mw"].dropna().to_numpy())
print(f"ANTES de quitar duplicados:   Mc={antes['mc']:.2f}  b={antes['b']:.3f}  "
      f"R2={antes['r2']:.3f}  N={antes['n_total']}")

df_sin_dup = df[~df["id"].isin(ids_usgs_duplicados)].copy()
print(f"\nCatalogo sin duplicados: N={len(df_sin_dup)}  "
      f"(se removieron {len(df) - len(df_sin_dup)} registros)")

despues = ajustar_gr(df_sin_dup["magnitude_mw"].dropna().to_numpy())
print(f"DESPUES de quitar duplicados: Mc={despues['mc']:.2f}  b={despues['b']:.3f}  "
      f"R2={despues['r2']:.3f}  N={despues['n_total']}")

print(f"\nDiferencia en b: {abs(antes['b'] - despues['b']):.4f}")

df_sin_dup.to_csv(
    DATOS / "SISMO NEIVA-CONSOLIDADO-260KM-HOMOGENEIZADO-MW-DEPURADO-GK-SINDUP.csv",
    index=False, encoding="utf-8-sig"
)
print(f"\n-> Catalogo sin duplicados guardado para referencia.")