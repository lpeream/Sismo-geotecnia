# -*- coding: utf-8 -*-
"""
Descarga los catalogos completos del SGC y del USGS alrededor del punto de
referencia (Parque Central Santander, Neiva) para los radios indicados.

Uso:  python descargar_catalogos.py 50 260
Salida: ../datos/datos originales/SISMO NEIVA-SGC-<R>KM-FULL.csv
        ../datos/datos originales/SISMO NEIVA-USGS-<R>KM.csv
"""

import csv
import json
import sys
import time
import urllib.request
from datetime import date

REF_LAT = 2.9262767
REF_LON = -75.2892111
SALIDA = "../datos/datos originales"

COLUMNAS_SGC = [
    "id", "status", "rms", "agency", "gap", "place", "local_time", "utc_time",
    "magnitude", "magnitude_error", "nst", "mag_type", "latitude",
    "latitude_error", "longitude", "longitude_error", "depth", "depth_error",
]


CABECERAS = {
    "User-Agent": (
        "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 "
        "(KHTML, like Gecko) Chrome/124.0 Safari/537.36"
    ),
    "Accept": "application/json, text/csv, */*",
    "Origin": "https://www.sgc.gov.co",
    "Referer": "https://www.sgc.gov.co/",
}


def pedir(url, intentos=4):
    for i in range(1, intentos + 1):
        try:
            req = urllib.request.Request(url, headers=CABECERAS)
            with urllib.request.urlopen(req, timeout=90) as r:
                return r.read()
        except Exception as e:
            if i == intentos:
                raise
            print(f"  reintento {i} tras error: {e}")
            time.sleep(3 * i)


def descargar_sgc(radio):
    url = (
        "https://apicatalogador.sgc.gov.co/api/events/search/"
        f"?lat_center={REF_LAT}&lon_center={REF_LON}&radius={radio}&page_size=1000"
    )
    eventos, total, pagina = [], 0, 1
    while url:
        datos = json.loads(pedir(url).decode("utf-8"))
        total = datos["count"]
        eventos.extend(datos["results"]["results"])
        print(f"  SGC {radio} km: pagina {pagina}, {len(eventos)}/{total}")
        url, pagina = datos["next"], pagina + 1
    assert len(eventos) == total, "El SGC devolvio menos eventos de los esperados"

    # La paginacion puede repetir una fila exacta si llega un sismo nuevo
    # mientras se descarga; solo se omite si es identica a la ya guardada.
    vistos, unicos = {}, []
    for e in eventos:
        if e["id"] in vistos and vistos[e["id"]] == e:
            print(f"  fila duplicada exacta omitida: {e['id']}")
            continue
        vistos[e["id"]] = e
        unicos.append(e)
    eventos = unicos

    archivo = f"{SALIDA}/SISMO NEIVA-SGC-{radio}KM-FULL.csv"
    with open(archivo, "w", encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(f, fieldnames=COLUMNAS_SGC, quoting=csv.QUOTE_MINIMAL)
        w.writeheader()
        for e in eventos:
            fila = {c: e.get(c) for c in COLUMNAS_SGC}
            for c in ("local_time", "utc_time"):
                if fila[c]:
                    fila[c] = fila[c].replace(" ", "T") + "+00:00"
            w.writerow(fila)
    print(f"  -> {archivo} ({len(eventos)} eventos)")


def descargar_usgs(radio):
    url = (
        "https://earthquake.usgs.gov/fdsnws/event/1/query?format=csv"
        f"&starttime=1900-01-01&endtime={date.today().isoformat()}"
        f"&latitude={REF_LAT}&longitude={REF_LON}&maxradiuskm={radio}"
        "&minmagnitude=0&orderby=time"
    )
    contenido = pedir(url)
    n = contenido.count(b"\n") - 1
    assert n < 20000, "USGS alcanzo el limite de 20000 eventos por consulta"
    archivo = f"{SALIDA}/SISMO NEIVA-USGS-{radio}KM.csv"
    with open(archivo, "wb") as f:
        f.write(contenido)
    print(f"  USGS {radio} km -> {archivo} ({n} eventos)")


if __name__ == "__main__":
    for r in sys.argv[1:]:
        descargar_sgc(int(r))
        descargar_usgs(int(r))
