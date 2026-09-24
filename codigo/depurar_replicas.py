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
UMBRAL_MW_TAREA = 3.0
DISTANCIA_REPLICA_TAREA_KM = 10.0
TIEMPO_REPLICA_TAREA_HORAS = 24.0

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
    columnas_requeridas = {"fecha_utc", "latitude", "longitude", "magnitude_mw"}
    faltantes = columnas_requeridas - set(df.columns)
    if faltantes:
        raise ValueError(f"Faltan columnas requeridas: {sorted(faltantes)}")

    resultado = df.copy()
    resultado["fecha_utc"] = pd.to_datetime(resultado["fecha_utc"], utc=True, errors="coerce")
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


def metricas_vecino_mas_cercano(df, alpha=1.0):
    """Calcula la métrica NN de Zaliapin et al. para cada evento y su vecino previo más cercano."""
    catalogo = df.copy()
    catalogo["fecha_utc"] = pd.to_datetime(catalogo["fecha_utc"], utc=True, errors="coerce")
    for columna in ("latitude", "longitude", "magnitude_mw"):
        catalogo[columna] = pd.to_numeric(catalogo[columna], errors="coerce")

    orden = catalogo.sort_values(["fecha_utc", "id"], na_position="last").index.tolist()
    metricas = []
    candidatos = []

    for indice_actual in orden:
        fecha_actual = catalogo.at[indice_actual, "fecha_utc"]
        lat_actual = catalogo.at[indice_actual, "latitude"]
        lon_actual = catalogo.at[indice_actual, "longitude"]
        if pd.isna(fecha_actual) or pd.isna(lat_actual) or pd.isna(lon_actual):
            continue

        prev_idx = []
        prev_dates = []
        prev_lats = []
        prev_lons = []

        for indice_prev in orden[: orden.index(indice_actual)]:
            fecha_prev = catalogo.at[indice_prev, "fecha_utc"]
            if pd.isna(fecha_prev):
                continue
            delta_t = (fecha_actual - fecha_prev).total_seconds() / 86400.0
            if delta_t <= 0:
                continue
            prev_idx.append(indice_prev)
            prev_dates.append(delta_t)
            prev_lats.append(catalogo.at[indice_prev, "latitude"])
            prev_lons.append(catalogo.at[indice_prev, "longitude"])

        if not prev_idx:
            continue

        distancias = distancia_km(
            lat_actual,
            lon_actual,
            np.asarray(prev_lats, dtype=float),
            np.asarray(prev_lons, dtype=float),
        )
        tiempos = np.asarray(prev_dates, dtype=float)
        metric = np.log10(np.maximum(distancias, 1.0)) + alpha * np.log10(np.maximum(tiempos, 1.0))
        mejor = int(np.argmin(metric))
        metricas.append(metric[mejor])
        candidatos.append(
            {
                "evento": indice_actual,
                "vecino": prev_idx[mejor],
                "delta_t_dias": tiempos[mejor],
                "distancia_km": float(distancias[mejor]),
                "metric_nn": float(metric[mejor]),
            }
        )

    if not candidatos:
        return pd.DataFrame(columns=["evento", "vecino", "delta_t_dias", "distancia_km", "metric_nn"]), None

    resumen = pd.DataFrame(candidatos)
    threshold = np.quantile(resumen["metric_nn"].dropna(), 0.10)
    return resumen, threshold


def declusterizar_nn(df, alpha=1.0):
    """Aplica la depuración por vecino más cercano siguiendo Zaliapin et al. (2008)."""
    columnas_requeridas = {"fecha_utc", "latitude", "longitude", "magnitude_mw"}
    faltantes = columnas_requeridas - set(df.columns)
    if faltantes:
        raise ValueError(f"Faltan columnas requeridas: {sorted(faltantes)}")

    resultado = df.copy()
    resultado["fecha_utc"] = pd.to_datetime(resultado["fecha_utc"], utc=True, errors="coerce")
    for columna in ("latitude", "longitude", "magnitude_mw"):
        resultado[columna] = pd.to_numeric(resultado[columna], errors="coerce")

    resultado["nn_tiempo_dias"] = pd.NA
    resultado["nn_distancia_km"] = pd.NA
    resultado["nn_metric"] = pd.NA
    resultado["nn_clasificacion"] = "principal"
    resultado["nn_mainshock_id"] = pd.NA
    resultado["nn_mainshock_mw"] = pd.NA

    nn_resumen, threshold = metricas_vecino_mas_cercano(resultado, alpha=alpha)
    if nn_resumen is None:
        resultado["fecha_utc"] = resultado["fecha_utc"].dt.strftime("%Y-%m-%dT%H:%M:%S.%fZ")
        return resultado.reset_index(drop=True), resultado

    threshold = threshold if threshold is not None else np.quantile(nn_resumen["metric_nn"], 0.10)
    excluidos = set()

    for _, fila in nn_resumen.iterrows():
        evento_idx = fila["evento"]
        vecino_idx = fila["vecino"]
        metric = fila["metric_nn"]
        if metric <= threshold:
            excluidos.add(evento_idx)
            resultado.at[evento_idx, "nn_clasificacion"] = "replica"
            resultado.at[evento_idx, "nn_tiempo_dias"] = fila["delta_t_dias"]
            resultado.at[evento_idx, "nn_distancia_km"] = fila["distancia_km"]
            resultado.at[evento_idx, "nn_metric"] = metric
            resultado.at[evento_idx, "nn_mainshock_id"] = resultado.at[vecino_idx, "id"]
            resultado.at[evento_idx, "nn_mainshock_mw"] = resultado.at[vecino_idx, "magnitude_mw"]

    resultado["fecha_utc"] = resultado["fecha_utc"].dt.strftime("%Y-%m-%dT%H:%M:%S.%fZ")
    return resultado.drop(index=sorted(excluidos)).reset_index(drop=True), resultado


def declusterizar_gkl(df):
    """Aplica la variante de ventanas enlazadas (GKL/GKL-b) usando la misma tabla GK."""
    columnas_requeridas = {"fecha_utc", "latitude", "longitude", "magnitude_mw"}
    faltantes = columnas_requeridas - set(df.columns)
    if faltantes:
        raise ValueError(f"Faltan columnas requeridas: {sorted(faltantes)}")

    resultado = df.copy()
    resultado["fecha_utc"] = pd.to_datetime(resultado["fecha_utc"], utc=True, errors="coerce")
    for columna in ("latitude", "longitude", "magnitude_mw"):
        resultado[columna] = pd.to_numeric(resultado[columna], errors="coerce")

    resultado["gkl_tiempo_dias"] = pd.NA
    resultado["gkl_distancia_km"] = pd.NA
    resultado["gkl_clasificacion"] = "principal"
    resultado["gkl_mainshock_id"] = pd.NA
    resultado["gkl_mainshock_mw"] = pd.NA

    orden = resultado.sort_values(["fecha_utc", "id"], na_position="last").index.tolist()
    excluidos = set()

    for indice_i in orden:
        if indice_i in excluidos:
            continue
        fecha_i = resultado.at[indice_i, "fecha_utc"]
        lat_i = resultado.at[indice_i, "latitude"]
        lon_i = resultado.at[indice_i, "longitude"]
        mw_i = resultado.at[indice_i, "magnitude_mw"]
        if pd.isna(fecha_i) or pd.isna(lat_i) or pd.isna(lon_i) or pd.isna(mw_i):
            continue

        padres = []
        for indice_j in orden[: orden.index(indice_i)]:
            if indice_j in excluidos:
                continue
            fecha_j = resultado.at[indice_j, "fecha_utc"]
            lat_j = resultado.at[indice_j, "latitude"]
            lon_j = resultado.at[indice_j, "longitude"]
            mw_j = resultado.at[indice_j, "magnitude_mw"]
            if pd.isna(fecha_j) or pd.isna(lat_j) or pd.isna(lon_j) or pd.isna(mw_j):
                continue
            delta_t = (fecha_i - fecha_j).total_seconds() / 86400.0
            if delta_t <= 0 or mw_j < mw_i:
                continue
            tiempo_j, distancia_j = ventanas_gardner_knopoff(mw_j)
            distancia = distancia_km(lat_i, lon_i, lat_j, lon_j)
            if delta_t <= tiempo_j and distancia <= distancia_j:
                padres.append((indice_j, delta_t, float(distancia), float(mw_j)))

        if not padres:
            continue

        indice_padre, delta_t, distancia, mw_padre = max(
            padres,
            key=lambda x: (x[3], -x[1]),
        )
        excluidos.add(indice_i)
        resultado.at[indice_i, "gkl_clasificacion"] = "replica"
        resultado.at[indice_i, "gkl_tiempo_dias"] = delta_t
        resultado.at[indice_i, "gkl_distancia_km"] = distancia
        resultado.at[indice_i, "gkl_mainshock_id"] = resultado.at[indice_padre, "id"]
        resultado.at[indice_i, "gkl_mainshock_mw"] = mw_padre

    resultado["fecha_utc"] = resultado["fecha_utc"].dt.strftime("%Y-%m-%dT%H:%M:%S.%fZ")
    return resultado.drop(index=sorted(excluidos)).reset_index(drop=True), resultado


def depurar_criterio_propio(
    df,
    umbral_mw=UMBRAL_MW_TAREA,
    distancia_maxima_km=DISTANCIA_REPLICA_TAREA_KM,
    tiempo_maximo_horas=TIEMPO_REPLICA_TAREA_HORAS,
):
    """Aplica el Criterio propio de depuracion."""
    columnas_requeridas = {"fecha_utc", "latitude", "longitude", "magnitude_mw"}
    faltantes = columnas_requeridas - set(df.columns)
    if faltantes:
        raise ValueError(f"Faltan columnas requeridas: {sorted(faltantes)}")

    resultado = df.copy()
    resultado["fecha_utc"] = pd.to_datetime(resultado["fecha_utc"], utc=True, errors="coerce")
    for columna in ("latitude", "longitude", "magnitude_mw"):
        resultado[columna] = pd.to_numeric(resultado[columna], errors="coerce")

    resultado["propio_clasificacion"] = "principal"
    resultado["propio_mainshock_id"] = pd.NA
    resultado["propio_mainshock_mw"] = pd.NA
    resultado["propio_distancia_km"] = pd.NA
    resultado["propio_motivo"] = pd.NA

    excluidos = set()
    magnitudes_bajas = resultado["magnitude_mw"].notna() & (
        resultado["magnitude_mw"] <= umbral_mw
    )
    for indice in resultado.index[magnitudes_bajas]:
        excluidos.add(indice)
        resultado.at[indice, "propio_clasificacion"] = "magnitud_baja"
        resultado.at[indice, "propio_motivo"] = f"Mw <= {umbral_mw:.1f}"

    elegibles = resultado.index.difference(excluidos)
    resultado["_dia_propio"] = resultado["fecha_utc"].dt.floor("D")
    grupos_por_dia = resultado.loc[elegibles].groupby("_dia_propio", sort=True).groups

    for indices_dia in grupos_por_dia.values():
        orden = resultado.loc[list(indices_dia)].sort_values(
            ["magnitude_mw", "id"], ascending=[False, True], na_position="last"
        ).index.tolist()
        for indice_principal in orden:
            if indice_principal in excluidos:
                continue
            fecha = resultado.at[indice_principal, "fecha_utc"]
            latitud = resultado.at[indice_principal, "latitude"]
            longitud = resultado.at[indice_principal, "longitude"]
            magnitud = resultado.at[indice_principal, "magnitude_mw"]
            if pd.isna(fecha) or pd.isna(latitud) or pd.isna(longitud) or pd.isna(magnitud):
                continue

            candidatos = [indice for indice in orden if indice != indice_principal and indice not in excluidos]
            if not candidatos:
                continue

            diferencias_horas = np.array(
                [
                    abs((resultado.at[indice, "fecha_utc"] - fecha).total_seconds()) / 3600.0
                    for indice in candidatos
                ]
            )
            candidatos = [
                indice for indice, diferencia in zip(candidatos, diferencias_horas)
                if diferencia <= tiempo_maximo_horas
                and resultado.at[indice, "magnitude_mw"] < magnitud
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
                if pd.isna(distancia) or distancia > distancia_maxima_km:
                    continue
                excluidos.add(indice_candidato)
                resultado.at[indice_candidato, "propio_clasificacion"] = "replica"
                resultado.at[indice_candidato, "propio_mainshock_id"] = resultado.at[
                    indice_principal, "id"
                ]
                resultado.at[indice_candidato, "propio_mainshock_mw"] = magnitud
                resultado.at[indice_candidato, "propio_distancia_km"] = float(distancia)
                resultado.at[indice_candidato, "propio_motivo"] = "mismo dia y cercania espacial"

    resultado = resultado.drop(columns="_dia_propio")
    resultado["fecha_utc"] = resultado["fecha_utc"].dt.strftime("%Y-%m-%dT%H:%M:%S.%fZ")
    return resultado.drop(index=sorted(excluidos)).reset_index(drop=True), resultado


def agregar_excluidos_por_magnitud(excluidos, catalogo, columna_clasificacion):
    """Agrega al reporte de exclusiones los eventos con Mw <= 3.0."""
    magnitudes = pd.to_numeric(catalogo["magnitude_mw"], errors="coerce")
    bajas = catalogo[magnitudes <= UMBRAL_MW_TAREA].copy()
    if bajas.empty:
        return excluidos
    bajas[columna_clasificacion] = "magnitud_baja"
    return pd.concat([bajas, excluidos], ignore_index=True, sort=False)


def procesar_radio(radio):
    entrada = DATOS / f"SISMO NEIVA-CONSOLIDADO-{radio}-HOMOGENEIZADO-MW.csv"

    catalogo_completo = pd.read_csv(entrada, encoding="utf-8-sig")
    magnitudes = pd.to_numeric(catalogo_completo["magnitude_mw"], errors="coerce")
    catalogo = catalogo_completo[magnitudes > UMBRAL_MW_TAREA].copy()

    depurado_gk, clasificado_gk = declusterizar(catalogo)
    depurado_gk.to_csv(
        DATOS / f"SISMO NEIVA-CONSOLIDADO-{radio}-HOMOGENEIZADO-MW-DEPURADO-GK.csv",
        index=False,
        encoding="utf-8-sig",
    )
    excluidos_gk = clasificado_gk[clasificado_gk["gk_clasificacion"] != "principal"]
    excluidos_gk = agregar_excluidos_por_magnitud(excluidos_gk, catalogo_completo, "gk_clasificacion")
    excluidos_gk.to_csv(
        DATOS / f"SISMO NEIVA-CONSOLIDADO-{radio}-REPLICAS-GK.csv",
        index=False,
        encoding="utf-8-sig",
    )

    depurado_nn, clasificado_nn = declusterizar_nn(catalogo)
    depurado_nn.to_csv(
        DATOS / f"SISMO NEIVA-CONSOLIDADO-{radio}-HOMOGENEIZADO-MW-DEPURADO-NN.csv",
        index=False,
        encoding="utf-8-sig",
    )
    excluidos_nn = clasificado_nn[clasificado_nn["nn_clasificacion"] != "principal"]
    excluidos_nn = agregar_excluidos_por_magnitud(excluidos_nn, catalogo_completo, "nn_clasificacion")
    excluidos_nn.to_csv(
        DATOS / f"SISMO NEIVA-CONSOLIDADO-{radio}-REPLICAS-NN.csv",
        index=False,
        encoding="utf-8-sig",
    )

    depurado_gkl, clasificado_gkl = declusterizar_gkl(catalogo)
    depurado_gkl.to_csv(
        DATOS / f"SISMO NEIVA-CONSOLIDADO-{radio}-HOMOGENEIZADO-MW-DEPURADO-GKL.csv",
        index=False,
        encoding="utf-8-sig",
    )
    excluidos_gkl = clasificado_gkl[clasificado_gkl["gkl_clasificacion"] != "principal"]
    excluidos_gkl = agregar_excluidos_por_magnitud(excluidos_gkl, catalogo_completo, "gkl_clasificacion")
    excluidos_gkl.to_csv(
        DATOS / f"SISMO NEIVA-CONSOLIDADO-{radio}-REPLICAS-GKL.csv",
        index=False,
        encoding="utf-8-sig",
    )

    depurado_tarea, clasificado_tarea = depurar_criterio_propio(catalogo)
    depurado_tarea.to_csv(
        DATOS / f"SISMO NEIVA-CONSOLIDADO-{radio}-HOMOGENEIZADO-MW-DEPURADO-PROPIO.csv",
        index=False,
        encoding="utf-8-sig",
    )
    excluidos_tarea = clasificado_tarea[clasificado_tarea["propio_clasificacion"] != "principal"]
    excluidos_tarea.to_csv(
        DATOS / f"SISMO NEIVA-CONSOLIDADO-{radio}-REPLICAS-PROPIO.csv",
        index=False,
        encoding="utf-8-sig",
    )

    print(
        f"{radio}: GK entrada={len(catalogo_completo)}  elegibles Mw>3={len(catalogo)}  "
        f"conservados={len(depurado_gk)}  "
        f"excluidos={len(excluidos_gk)}"
    )
    print(excluidos_gk["gk_clasificacion"].value_counts().to_string())
    print(
        f"{radio}: NN entrada={len(catalogo_completo)}  elegibles Mw>3={len(catalogo)}  "
        f"conservados={len(depurado_nn)}  "
        f"excluidos={len(excluidos_nn)}"
    )
    print(excluidos_nn["nn_clasificacion"].value_counts().to_string())
    print(
        f"{radio}: GKL entrada={len(catalogo_completo)}  elegibles Mw>3={len(catalogo)}  "
        f"conservados={len(depurado_gkl)}  "
        f"excluidos={len(excluidos_gkl)}"
    )
    print(excluidos_gkl["gkl_clasificacion"].value_counts().to_string())
    print(
        f"{radio}: PROPIO entrada={len(catalogo_completo)}  elegibles Mw>3={len(catalogo)}  "
        f"conservados={len(depurado_tarea)}  "
        f"excluidos={len(excluidos_tarea)}"
    )
    print(excluidos_tarea["propio_clasificacion"].value_counts().to_string())


for radio in ("50KM", "260KM"):
    procesar_radio(radio)