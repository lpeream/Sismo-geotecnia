# -*- coding: utf-8 -*-
"""
Estima la magnitud de completitud (Mc) por separado para cada fuente
(SGC, USGS) y por periodos de tiempo dentro del SGC, para ambos radios
(50km y 260km), usando el metodo de maxima curvatura.
Corre sobre el catalogo YA depurado (post Gardner-Knopoff, con el bug
de fechas ya corregido).
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DATOS = ROOT / "datos" / "datos depurados"
FIGURAS = ROOT / "figuras"
FIGURAS.mkdir(exist_ok=True)
ANCHO_BIN = 0.1


def ajustar_gr(magnitudes):
    if len(magnitudes) < 10:
        return None

    minimo = np.floor(magnitudes.min() / ANCHO_BIN) * ANCHO_BIN
    maximo = np.ceil(magnitudes.max() / ANCHO_BIN) * ANCHO_BIN + ANCHO_BIN
    bordes = np.arange(minimo, maximo + ANCHO_BIN / 2, ANCHO_BIN)
    frecuencias, _ = np.histogram(magnitudes, bins=bordes)
    centros = bordes[:-1] + ANCHO_BIN / 2
    indice_mc = int(np.argmax(frecuencias))
    mc = float(centros[indice_mc])

    valores = np.sort(magnitudes[magnitudes >= mc])
    magnitudes_unicas = np.unique(valores)
    if len(magnitudes_unicas) < 3:
        return {"mc": mc, "b": np.nan, "r2": np.nan, "n_total": len(magnitudes), "n_ajuste": len(valores)}

    acumuladas = np.array([(valores >= v).sum() for v in magnitudes_unicas])
    mascara = acumuladas >= 2
    x = magnitudes_unicas[mascara]
    y = np.log10(acumuladas[mascara])
    if len(x) < 2:
        return {"mc": mc, "b": np.nan, "r2": np.nan, "n_total": len(magnitudes), "n_ajuste": len(valores)}

    pendiente, intercepto = np.polyfit(x, y, 1)
    prediccion = intercepto + pendiente * x
    ss_res = np.sum((y - prediccion) ** 2)
    ss_tot = np.sum((y - y.mean()) ** 2)
    r2 = 1 - ss_res / ss_tot if ss_tot else np.nan

    return {
        "mc": mc, "a": intercepto, "b": -pendiente, "r2": r2,
        "n_total": len(magnitudes), "n_ajuste": len(valores),
    }


resultados = []

for radio in ["50KM", "260KM"]:
    archivo = DATOS / f"SISMO NEIVA-CONSOLIDADO-{radio}-HOMOGENEIZADO-MW-DEPURADO-GK.csv"
    df = pd.read_csv(archivo, encoding="utf-8-sig")
    df["magnitude_mw"] = pd.to_numeric(df["magnitude_mw"], errors="coerce")
    df["fecha_utc"] = pd.to_datetime(df["fecha_utc"], errors="coerce", format="ISO8601")
    df["anio"] = df["fecha_utc"].dt.year

    print("=" * 60)
    print(f"RADIO: {radio}  (fechas no parseadas: {df['fecha_utc'].isna().sum()} de {len(df)})")
    print("=" * 60)

    print(f"\n--- Completitud por fuente ({radio}) ---")
    for fuente in ["SGC", "USGS"]:
        subset = df[df["fuente"] == fuente].dropna(subset=["magnitude_mw"])
        ajuste = ajustar_gr(subset["magnitude_mw"].to_numpy())
        if ajuste:
            resultados.append({"grupo": f"{radio}_{fuente}", **ajuste})
            print(f"{fuente}: Mc={ajuste['mc']:.2f}  b={ajuste.get('b', float('nan')):.3f}  "
                  f"R2={ajuste.get('r2', float('nan')):.3f}  N={ajuste['n_total']}")
        else:
            print(f"{fuente}: datos insuficientes (N={len(subset)})")

    print(f"\n--- Completitud por periodo, solo SGC ({radio}) ---")
    sgc = df[df["fuente"] == "SGC"].dropna(subset=["magnitude_mw", "anio"])
    if len(sgc) > 0:
        anio_min, anio_max = int(sgc["anio"].min()), int(sgc["anio"].max())
        print(f"Rango de anios en SGC: {anio_min} - {anio_max}")
        cortes = np.linspace(anio_min, anio_max + 1, 4)
        for i in range(3):
            ini, fin = int(cortes[i]), int(cortes[i + 1])
            subset = sgc[(sgc["anio"] >= ini) & (sgc["anio"] < fin)]
            ajuste = ajustar_gr(subset["magnitude_mw"].to_numpy())
            etiqueta = f"{radio}_SGC_{ini}-{fin-1}"
            if ajuste:
                resultados.append({"grupo": etiqueta, **ajuste})
                print(f"{etiqueta}: Mc={ajuste['mc']:.2f}  b={ajuste.get('b', float('nan')):.3f}  "
                      f"R2={ajuste.get('r2', float('nan')):.3f}  N={ajuste['n_total']}")
            else:
                print(f"{etiqueta}: datos insuficientes")
    print()

tabla = pd.DataFrame(resultados)
salida = DATOS / "completitud_por_fuente_periodo.csv"
tabla.to_csv(salida, index=False, encoding="utf-8-sig")
print(f"-> Guardado en: {salida}")

fig, ax = plt.subplots(figsize=(11, 5.5))
ax.bar(tabla["grupo"], tabla["mc"], color="#1f5d75")
ax.set_ylabel(r"Magnitud de completitud $M_c$")
ax.set_title("Completitud por radio, fuente y periodo (SGC)")
ax.tick_params(axis="x", rotation=45, labelsize=8)
fig.tight_layout()
fig.savefig(FIGURAS / "completitud_por_fuente.png", dpi=220, bbox_inches="tight")
print("-> Figura guardada en figuras/completitud_por_fuente.png")