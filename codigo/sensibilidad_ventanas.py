# -*- coding: utf-8 -*-
"""
Compara la depuracion Gardner-Knopoff contra la ventana alternativa de
Uhrhammer (1986), usando el catalogo YA homogeneizado (antes de GK), para
evaluar la sensibilidad del numero de exclusiones a la tabla de ventanas
elegida.
"""

from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DATOS = ROOT / "datos" / "datos depurados"


def ventanas_uhrhammer(magnitud):
    """Uhrhammer (1986): formulas continuas, no tabla interpolada."""
    tiempo = np.exp(-2.87 + 0.235 * magnitud)
    distancia = np.exp(-1.024 + 0.804 * magnitud)
    return tiempo, distancia


def distancia_km(latitud_1, longitud_1, latitud_2, longitud_2):
    radio_tierra = 6371.0
    latitud_1 = np.radians(latitud_1)
    latitud_2 = np.radians(latitud_2)
    diferencia_latitud = latitud_2 - latitud_1
    diferencia_longitud = np.radians(longitud_2 - longitud_1)
    termino = (
        np.sin(diferencia_latitud / 2) ** 2
        + np.cos(latitud_1) * np.cos(latitud_2) * np.sin(diferencia_longitud / 2) ** 2
    )
    return 2 * radio_tierra * np.arcsin(np.sqrt(termino))


def declusterizar_uhrhammer(df):
    resultado = df.copy()
    resultado["fecha_utc"] = pd.to_datetime(resultado["fecha_utc"], utc=True, errors="coerce", format="ISO8601")
    for columna in ("latitude", "longitude", "magnitude_mw"):
        resultado[columna] = pd.to_numeric(resultado[columna], errors="coerce")

    resultado["gk_clasificacion"] = "principal"

    orden = resultado.sort_values(
        ["magnitude_mw", "fecha_utc", "id"], ascending=[False, True, True], na_position="last"
    ).index.tolist()
    excluidos = set()

    for indice_principal in orden:
        if indice_principal in excluidos:
            continue
        magnitud = resultado.at[indice_principal, "magnitude_mw"]
        fecha = resultado.at[indice_principal, "fecha_utc"]
        latitud = resultado.at[indice_principal, "latitude"]
        longitud = resultado.at[indice_principal, "longitude"]
        if pd.isna(magnitud) or pd.isna(fecha) or pd.isna(latitud) or pd.isna(longitud):
            continue

        tiempo_dias, distancia_maxima = ventanas_uhrhammer(magnitud)

        fechas = resultado["fecha_utc"]
        candidatos = resultado.index[
            fechas.between(
                fecha - pd.Timedelta(days=tiempo_dias),
                fecha + pd.Timedelta(days=tiempo_dias),
                inclusive="both",
            )
        ]
        candidatos = [i for i in candidatos if i != indice_principal and i not in excluidos]
        if not candidatos:
            continue

        distancias = distancia_km(
            latitud, longitud,
            resultado.loc[candidatos, "latitude"].to_numpy(),
            resultado.loc[candidatos, "longitude"].to_numpy(),
        )
        for indice_candidato, distancia in zip(candidatos, distancias):
            if pd.isna(distancia) or distancia > distancia_maxima:
                continue
            fecha_candidato = resultado.at[indice_candidato, "fecha_utc"]
            excluidos.add(indice_candidato)
            resultado.at[indice_candidato, "gk_clasificacion"] = (
                "replica" if fecha_candidato >= fecha else "precursor"
            )

    return resultado.drop(index=sorted(excluidos)).reset_index(drop=True), resultado, excluidos


resumen = []
for radio in ("50KM", "260KM"):
    entrada = DATOS / f"SISMO NEIVA-CONSOLIDADO-{radio}-HOMOGENEIZADO-MW.csv"
    catalogo = pd.read_csv(entrada, encoding="utf-8-sig")
    depurado, clasificado, excluidos = declusterizar_uhrhammer(catalogo)
    depurado.to_csv(DATOS / f"SISMO NEIVA-CONSOLIDADO-{radio}-UHRHAMMER-DEPURADO.csv", index=False, encoding="utf-8-sig")
    n_excluidos = len(excluidos)
    n_replicas = (clasificado.loc[list(excluidos), "gk_clasificacion"] == "replica").sum()
    n_precursores = (clasificado.loc[list(excluidos), "gk_clasificacion"] == "precursor").sum()
    print(f"{radio} (Uhrhammer): entrada={len(catalogo)}  conservados={len(depurado)}  "
          f"excluidos={n_excluidos}  (replicas={n_replicas}, precursores={n_precursores})")
    resumen.append({"radio": radio, "metodo": "Uhrhammer", "entrada": len(catalogo),
                     "conservados": len(depurado), "excluidos": n_excluidos,
                     "replicas": n_replicas, "precursores": n_precursores})

tabla = pd.DataFrame(resumen)
tabla.to_csv(DATOS / "sensibilidad_ventanas_uhrhammer.csv", index=False, encoding="utf-8-sig")
print(f"\n-> Guardado en: {DATOS / 'sensibilidad_ventanas_uhrhammer.csv'}")