"""Elimina replicas y precursores mediante ventanas Gardner-Knopoff.

Los eventos se procesan de mayor a menor magnitud. Cada evento que no haya
sido asignado y que caiga dentro de la ventana espacio-tiempo del evento
principal se marca como replica o precursor y se excluye del catalogo.
"""

from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DATOS = ROOT / "datos" / "datos depurados"

# Tabla de ventanas Gardner-Knopoff (tiempo en dias, distancia en km).
MAGNITUDES_GK = np.array(
    [2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0]
)
TIEMPO_GK = np.array(
    [1, 1, 2, 3, 4, 6, 10, 15, 22, 32, 42, 83, 155, 290, 510], dtype=float
)
DISTANCIA_GK = np.array(
    [5, 5, 5, 5, 10, 15, 20, 30, 50, 70, 90, 120, 150, 180, 210], dtype=float
)


def ventanas_gardner_knopoff(magnitud):
    """Interpola las ventanas GK y limita los extremos a la tabla publicada."""
    magnitud_limitada = np.clip(magnitud, MAGNITUDES_GK[0], MAGNITUDES_GK[-1])
    tiempo = np.interp(magnitud_limitada, MAGNITUDES_GK, TIEMPO_GK)
    distancia = np.interp(magnitud_limitada, MAGNITUDES_GK, DISTANCIA_GK)
    return tiempo, distancia


def distancia_km(latitud_1, longitud_1, latitud_2, longitud_2):
    """Calcula la distancia de gran circulo con la formula de Haversine."""
    radio_tierra = 6371.0
    latitud_1 = np.radians(latitud_1)
    latitud_2 = np.radians(latitud_2)
    diferencia_latitud = latitud_2 - latitud_1
    diferencia_longitud = np.radians(longitud_2 - longitud_1)
    termino = (
        np.sin(diferencia_latitud / 2) ** 2
        + np.cos(latitud_1)
        * np.cos(latitud_2)
        * np.sin(diferencia_longitud / 2) ** 2
    )
    return 2 * radio_tierra * np.arcsin(np.sqrt(termino))


def declusterizar(df):
    """Devuelve el catalogo depurado y las ventanas aplicadas a cada evento."""
    columnas_requeridas = {
        "fecha_utc", "latitude", "longitude", "magnitude_mw"
    }
    faltantes = columnas_requeridas - set(df.columns)
    if faltantes:
        raise ValueError(f"Faltan columnas requeridas: {sorted(faltantes)}")

    resultado = df.copy()
    resultado["fecha_utc"] = pd.to_datetime(resultado["fecha_utc"], utc=True, errors="coerce", format="ISO8601")
    for columna in ("latitude", "longitude", "magnitude_mw"):
        resultado[columna] = pd.to_numeric(resultado[columna], errors="coerce")

    resultado["gk_tiempo_dias"] = pd.NA
    resultado["gk_distancia_km"] = pd.NA
    resultado["gk_clasificacion"] = "principal"
    resultado["gk_mainshock_id"] = pd.NA
    resultado["gk_mainshock_mw"] = pd.NA

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

        tiempo_dias, distancia_maxima = ventanas_gardner_knopoff(magnitud)
        resultado.at[indice_principal, "gk_tiempo_dias"] = tiempo_dias
        resultado.at[indice_principal, "gk_distancia_km"] = distancia_maxima

        fechas = resultado["fecha_utc"]
        candidatos = resultado.index[
            fechas.between(
                fecha - pd.Timedelta(days=tiempo_dias),
                fecha + pd.Timedelta(days=tiempo_dias),
                inclusive="both",
            )
        ]
        candidatos = [
            indice for indice in candidatos
            if indice != indice_principal and indice not in excluidos
        ]
        if not candidatos:
            continue

        distancias = distancia_km(
            latitud,
            longitud,
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
            resultado.at[indice_candidato, "gk_mainshock_id"] = resultado.at[
                indice_principal, "id"
            ]
            resultado.at[indice_candidato, "gk_mainshock_mw"] = magnitud


    resultado["fecha_utc"] = resultado["fecha_utc"].dt.strftime("%Y-%m-%dT%H:%M:%S.%fZ")
    return resultado.drop(index=sorted(excluidos)).reset_index(drop=True), resultado


def procesar_radio(radio):
    entrada = DATOS / f"SISMO NEIVA-CONSOLIDADO-{radio}-HOMOGENEIZADO-MW.csv"
    salida = DATOS / f"SISMO NEIVA-CONSOLIDADO-{radio}-HOMOGENEIZADO-MW-DEPURADO-GK.csv"
    trazabilidad = DATOS / f"SISMO NEIVA-CONSOLIDADO-{radio}-REPLICAS-GK.csv"

    catalogo = pd.read_csv(entrada, encoding="utf-8-sig")
    depurado, clasificado = declusterizar(catalogo)
    depurado.to_csv(salida, index=False, encoding="utf-8-sig")
    excluidos = clasificado[clasificado["gk_clasificacion"] != "principal"]
    excluidos.to_csv(trazabilidad, index=False, encoding="utf-8-sig")


    print(
        f"{radio}: entrada={len(catalogo)}  conservados={len(depurado)}  "
        f"excluidos={len(excluidos)}"
    )
    print(excluidos["gk_clasificacion"].value_counts().to_string())


for radio in ("50KM", "260KM"):
    procesar_radio(radio)