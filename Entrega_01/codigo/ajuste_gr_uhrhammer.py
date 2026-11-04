# -*- coding: utf-8 -*-
"""
Ajusta Gutenberg-Richter sobre el catalogo depurado con las ventanas de
Uhrhammer (1986), para comparar el parametro b contra el obtenido con
Gardner-Knopoff.
"""

from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DATOS = ROOT / "datos" / "datos depurados"
ANCHO_BIN = 0.1


def ajustar_gr(magnitudes):
    if len(magnitudes) < 10:
        return None
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
    ss_res = np.sum((y - prediccion) ** 2)
    ss_tot = np.sum((y - y.mean()) ** 2)
    r2 = 1 - ss_res / ss_tot if ss_tot else np.nan

    return {"mc": mc, "b": -pendiente, "r2": r2, "n_total": len(magnitudes), "n_ajuste": len(valores)}


print("=" * 60)
print("COMPARACION: Gardner-Knopoff vs Uhrhammer")
print("=" * 60)

for radio in ["50KM", "260KM"]:
    archivo = DATOS / f"SISMO NEIVA-CONSOLIDADO-{radio}-UHRHAMMER-DEPURADO.csv"
    df = pd.read_csv(archivo, encoding="utf-8-sig")
    df["magnitude_mw"] = pd.to_numeric(df["magnitude_mw"], errors="coerce")
    magnitudes = df["magnitude_mw"].dropna().to_numpy()
    ajuste = ajustar_gr(magnitudes)
    print(f"\n{radio} (Uhrhammer): Mc={ajuste['mc']:.2f}  b={ajuste['b']:.3f}  "
          f"R2={ajuste['r2']:.3f}  N_total={ajuste['n_total']}  N_ajuste={ajuste['n_ajuste']}")