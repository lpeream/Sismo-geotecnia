"""Grafica y ajusta la distribucion frecuencia-magnitud de los catalogos Mw.

El ajuste sigue la relacion Gutenberg-Richter:
    log10 N(Mw >= M) = a - b M

La magnitud de completitud Mc se estima por maxima curvatura del histograma
con bins de 0.1 unidades. Los resultados son exploratorios porque los
catalogos combinan fuentes, escalas y conversiones con distinta incertidumbre.

Se repite el analisis para los catalogos depurados con Gardner-Knopoff,
ventanas enlazadas, vecino mas cercano y el Criterio propio.
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
METODOS = ("GK", "GKL", "NN", "PROPIO")
NOMBRES_METODOS = {
    "GK": "Gardner-Knopoff",
    "GKL": "Ventanas enlazadas (GKL)",
    "NN": "Vecino más cercano",
    "PROPIO": "Criterio propio",
}


def cargar_magnitudes(radio, metodo):
    archivo = DATOS / f"SISMO NEIVA-CONSOLIDADO-{radio}-HOMOGENEIZADO-MW-DEPURADO-{metodo}.csv"
    df = pd.read_csv(archivo, encoding="utf-8-sig")
    magnitudes = pd.to_numeric(df["magnitude_mw"], errors="coerce").dropna()
    return df, magnitudes.to_numpy()


def ajustar_gr(magnitudes):
    minimo = np.floor(magnitudes.min() / ANCHO_BIN) * ANCHO_BIN
    maximo = np.ceil(magnitudes.max() / ANCHO_BIN) * ANCHO_BIN + ANCHO_BIN
    bordes = np.arange(minimo, maximo + ANCHO_BIN / 2, ANCHO_BIN)
    frecuencias, _ = np.histogram(magnitudes, bins=bordes)
    centros = bordes[:-1] + ANCHO_BIN / 2
    indice_mc = int(np.argmax(frecuencias))
    mc = float(centros[indice_mc])
    umbral = mc

    valores = np.sort(magnitudes[magnitudes >= umbral])
    magnitudes_unicas = np.unique(valores)
    acumuladas = np.array([(valores >= valor).sum() for valor in magnitudes_unicas])
    mascara = acumuladas >= 2
    x = magnitudes_unicas[mascara]
    y = np.log10(acumuladas[mascara])
    pendiente, intercepto = np.polyfit(x, y, 1)
    prediccion = intercepto + pendiente * x
    ss_res = np.sum((y - prediccion) ** 2)
    ss_tot = np.sum((y - y.mean()) ** 2)
    r2 = 1 - ss_res / ss_tot if ss_tot else np.nan
    b = -pendiente
    a = intercepto

    return {
        "mc": mc,
        "umbral": umbral,
        "a": a,
        "b": b,
        "r2": r2,
        "n_total": len(magnitudes),
        "n_ajuste": len(valores),
        "centros": centros,
        "frecuencias": frecuencias,
        "x": x,
        "acumuladas": acumuladas[mascara],
    }


def graficar(radio, metodo):
    df, magnitudes = cargar_magnitudes(radio, metodo)
    archivo_replicas = DATOS / f"SISMO NEIVA-CONSOLIDADO-{radio}-REPLICAS-{metodo}.csv"
    n_replicas = len(pd.read_csv(archivo_replicas, encoding="utf-8-sig"))
    ajuste = ajustar_gr(magnitudes)
    nombre_metodo = NOMBRES_METODOS.get(metodo, metodo)
    figura, ejes = plt.subplots(1, 2, figsize=(13, 5.5))
    color = "#1f5d75"
    color_ajuste = "#b24732"

    ejes[0].bar(
        ajuste["centros"],
        ajuste["frecuencias"],
        width=ANCHO_BIN * 0.9,
        color=color,
        edgecolor="white",
        linewidth=0.25,
    )
    ejes[0].axvline(
        ajuste["mc"],
        color=color_ajuste,
        linestyle="--",
        linewidth=1.5,
        label=fr"$M_c={ajuste['mc']:.1f}$",
    )
    ejes[0].set_title(f"Frecuencia por magnitud - {radio.replace('KM', ' km')} ({nombre_metodo})")
    ejes[0].set_xlabel(r"Magnitud homogenizada $M_w$")
    ejes[0].set_ylabel("Número de eventos")
    ejes[0].legend(frameon=False)
    ejes[0].grid(axis="y", alpha=0.25)

    ejes[1].scatter(
        ajuste["x"],
        np.log10(ajuste["acumuladas"]),
        s=18,
        color=color,
        alpha=0.8,
        label="Catalogo acumulado",
    )
    limite = np.array([ajuste["umbral"], magnitudes.max()])
    ejes[1].plot(
        limite,
        ajuste["a"] - ajuste["b"] * limite,
        color=color_ajuste,
        linewidth=2,
        label=fr"$\log_{{10}}N={ajuste['a']:.2f}-{ajuste['b']:.2f}M_w$",
    )
    ejes[1].axvline(ajuste["mc"], color=color_ajuste, linestyle="--", linewidth=1)
    ejes[1].set_title(f"Ajuste Gutenberg-Richter ({nombre_metodo})")
    ejes[1].set_xlabel(r"Magnitud $M_w$")
    ejes[1].set_ylabel(r"$\log_{10} N(M_w \geq M)$")
    ejes[1].legend(frameon=False, loc="upper right")
    ejes[1].grid(alpha=0.25)

    estado = df["conversion_status"].value_counts().to_dict()
    texto = (
        f"N valido = {ajuste['n_total']:,}\n"
        f"N ajuste = {ajuste['n_ajuste']:,}\n"
        f"b = {ajuste['b']:.3f}; R2 = {ajuste['r2']:.3f}\n"
        f"Extrapoladas = {estado.get('extrapolada', 0):,}\n"
        f"Replicas/ precursores excluidos = {n_replicas:,}"
    )
    figura.text(0.5, 0.01, texto, ha="center", va="bottom", fontsize=9)
    figura.tight_layout(rect=(0, 0.05, 1, 1))
    salida = FIGURAS / f"frecuencia_magnitud_{radio.lower()}_{metodo.lower()}.png"
    figura.savefig(salida, dpi=220, bbox_inches="tight")
    plt.close(figura)
    return ajuste


resultados = []
for metodo in METODOS:
    for radio in ("50KM", "260KM"):
        ajuste = graficar(radio, metodo)
        resultados.append(
            {
                "metodo": metodo,
                "radio": radio,
                "Mc": ajuste["mc"],
                "a": ajuste["a"],
                "b": ajuste["b"],
                "R2": ajuste["r2"],
                "N_total": ajuste["n_total"],
                "N_ajuste": ajuste["n_ajuste"],
            }
        )

pd.DataFrame(resultados).to_csv(
    DATOS / "ajuste_gutenberg_richter.csv", index=False, encoding="utf-8-sig"
)
for resultado in resultados:
    print(resultado)
