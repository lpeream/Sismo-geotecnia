"""PDF de lectura de la onda S, canal por canal (figuras/lectura_downhole_S/).

Incluye instrucciones, un ejemplo resuelto, las decisiones del procesamiento,
una página por canal (golpes X y Y con la llegada leída) y la gráfica tz vs z.
Los tiempos se toman de procesamiento_downhole.py (T_S_MS).
"""
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import MultipleLocator

ROOT = Path(__file__).resolve().parents[1]   # carpeta Entrega_03
sys.path.insert(0, str(ROOT / "codigo"))
from procesamiento_downhole import T_S_MS, read_records  # noqa: E402

out = ROOT / "figuras" / "lectura_downhole_S"
out.mkdir(parents=True, exist_ok=True)
_t, _p, SX, SY = read_records()
t = _t * 1000
A4 = (11.69, 8.27)
T = np.asarray(T_S_MS, dtype=float)[:30]
Z = np.arange(1, 31.); R_ = np.hypot(Z, 3); TZ = T * Z / R_
pd.DataFrame({"canal": Z.astype(int), "z_m": Z, "t_S_ms": T, "R_m": R_.round(3), "tz_S_ms": TZ.round(2),
              "lectura": ["visual: inicio de la S fuerte"] * 17 + ["relativa al canal 17 (correlación)"] * 13}
             ).to_csv(ROOT / "datos" / "resultados" / "downhole_tiempos_S_lectura.csv", index=False,
                      encoding="utf-8-sig")


def text_page(pdf, title, blocks):
    fig = plt.figure(figsize=A4); fig.text(0.06, 0.93, title, fontsize=18, weight="bold")
    y = 0.86
    for head, lines in blocks:
        if head: fig.text(0.06, y, head, fontsize=13, weight="bold", color="#1f4e79"); y -= 0.045
        for ln in lines:
            fig.text(0.08, y, ln, fontsize=11.5); y -= 0.037
        y -= 0.02
    pdf.savefig(fig); plt.close(fig)

def channel_ax(ax, k, x0, x1):
    m = (t >= x0) & (t <= x1)
    ax.plot(t[m], SX[k][m], color="#1f4e79", lw=1.3, label="Golpe X")
    ax.plot(t[m], SY[k][m], color="#c0392b", lw=1.3, label="Golpe Y")
    ax.axhline(0, color="k", lw=.5); ax.set_xlim(x0, x1)
    ax.xaxis.set_major_locator(MultipleLocator(10)); ax.xaxis.set_minor_locator(MultipleLocator(1))
    ax.grid(which="major", alpha=.6); ax.grid(which="minor", alpha=.2)
    ax.set_xlabel("Tiempo (ms)"); ax.set_ylabel("Amplitud"); ax.legend(loc="upper left")

with PdfPages(out / "canales_1_a_30.pdf") as pdf:
    text_page(pdf, "Lectura de la onda S en el down hole - instrucciones", [
        ("Qué hay en el Excel (hoja VS)", [
            "Cada canal (1 m de profundidad) tiene dos golpes: X y Y. Cada golpe tiene dos columnas:",
            "  - Tiempo: el reloj del registro (0; 0,5; 1,0 ms ...). Es igual para todos; por sí solo no dice cuándo llegó la onda.",
            "  - Amplitud: cuánto se movió el geófono en cada instante.",
            "Para saber cuándo llegó la onda se grafica amplitud contra tiempo (las páginas siguientes)."]),
        ("Qué se ve en cada gráfica", [
            "1. Al principio la señal está quieta (cerca de cero).",
            "2. Luego aparecen ondulaciones pequeñas: es la onda P (llega primero, viaja más rápido). NO se lee.",
            "3. Después X y Y se separan en espejo (una sube y la otra baja) y con amplitud grande: es la onda S.",
            "   En los canales de arriba la señal se recorta (queda plana arriba y abajo): lo que importa es dónde empieza."]),
        ("Por qué X y Y", [
            "X y Y son golpes horizontales en sentidos contrarios. La onda S cambia de sentido con el golpe,",
            "por eso aparece en espejo. Así se distingue la S de la P y del ruido."]),
        ("Qué se anota", [
            "Solo el TIEMPO (eje horizontal) donde empieza la S, en milisegundos. La amplitud no entra en las fórmulas;",
            "sirve solo para ver en la gráfica cuándo llegó la onda. Si X y Y difieren 1-2 ms, tome el primero o el promedio",
            "y use el mismo criterio en todos los canales. El tiempo debe aumentar con la profundidad."]),
    ])
    # página de ejemplo con el canal 5 anotado
    fig, ax = plt.subplots(figsize=A4)
    channel_ax(ax, 4, 0, 80)
    ax.axvspan(9, 37, color="#c08a2b", alpha=.15)
    ax.text(23, 0.0045, "Onda P (pequeña)\nNO se lee", ha="center", fontsize=13, color="#8a5a00")
    ax.axvline(41, color="green", lw=2.5)
    ax.annotate("Llegada de la onda S\nt = 41 ms\n(X y Y se separan en espejo\ncon amplitud grande)", xy=(41, -0.002),
                xytext=(52, -0.0045), fontsize=13, color="green", arrowprops=dict(arrowstyle="->", color="green", lw=2))
    ax.set_title("Ejemplo: canal 5 (z = 5 m)", fontsize=15)
    fig.tight_layout(rect=(0, 0.12, 1, 1))
    fig.text(0.06, 0.06, "Con t = 41 ms:  R = √(5² + 3²) = 5,83 m   →   tz = t·z/R = 41·5/5,83 = 35,2 ms", fontsize=13)
    pdf.savefig(fig); plt.close(fig)
    text_page(pdf, "Del tiempo leído a Vs", [
        ("Paso 1. Leer t en cada canal (páginas siguientes)", ["Anotarlo en la columna amarilla 't S' de calculo_vs_downhole.xlsx (hoja CALCULO_VS)."]),
        ("Paso 2. Corregir a tiempo vertical", [
            "La fuente está a 3 m del pozo; la onda viaja en diagonal. R = √(z² + 3²)      tz = t · z / R",
            "La corrección pesa en los primeros metros (canal 2: tz ≈ 0,55·t) y casi nada abajo (canal 30: tz ≈ 0,995·t)."]),
        ("Paso 3. Graficar tz contra z", [
            "Los puntos forman tramos rectos; cada tramo es una capa. Donde cambia la inclinación hay cambio de capa.",
            "Tramo donde tz crece mucho por metro = suelo blando (Vs baja). Tramo donde crece poco = suelo rígido (Vs alta)."]),
        ("Paso 4. Velocidad de cada capa", [
            "Vs = 1 / pendiente de la recta tz contra z = Δz / Δtz.   En Excel: =1000/PENDIENTE(tz; z)  (tz en ms).",
            "Se usa la recta con varios canales, no metro a metro: un error de 1-2 ms en un canal cambia mucho Vs entre",
            "dos canales vecinos, pero se compensa en la recta."]),
        ("Paso 5. Resultados del perfil", [
            "Vs30 = 30 / Σ(h / Vs)          Gmax = (γ / g) · Vs²",
            "En la hoja CALCULO_VS esto se calcula solo al llenar los tiempos y los canales de cada capa."]),
    ])
    text_page(pdf, "Decisiones del procesamiento", [
        ("Lecturas", [
            "Canales 1 a 17: lectura visual en el inicio de la onda S fuerte (línea verde en cada página).",
            "Canales 18 a 30: el inicio no es claro; t = t(canal 17) + retardos de la correlación cruzada entre canales",
            "vecinos (marcados como 'lectura relativa')."]),
        ("Canales que no se usan", [
            "Canales 1 y 2 (1 y 2 m): campo cercano. Están a menos de una longitud de onda de la fuente",
            "(λ = Vs/f ≈ 145/50 ≈ 3 m) y la señal está saturada. Se muestran pero no entran en las rectas."]),
        ("Capas (rectas de tz contra z)", [
            "Arena: 0-6 m (canales 3-6), 6-11 m (6-11) y 11-16,5 m (11-16). Arcilla: una recta, 16,5-30 m (17-30).",
            "El cambio de pendiente está entre los canales 16 y 17: el límite arena-arcilla se toma en 16,5 m",
            "(referencia: 16,4 m) y ninguna recta lo cruza."]),
        ("Salto entre 17 y 20 m", [
            "Del canal 17 al 18 la correlación daba un salto de ~10 ms (un ciclo de la onda); el canal 18 se ajustó a 127 ms.",
            "en la correlación, no un cambio de suelo (la muestra M22 y la columna resonante de 20 m siguen siendo la",
            "misma arcilla). La recta de 17 a 30 m reparte el salto."]),
        ("Resultado", [
            "Arena 145 / 148 / 160 m/s (151 m/s el estrato) y arcilla 380 m/s.   Vs30 = 207 m/s: perfil tipo D (NSR-10)."]),
    ])
    for k in range(30):
        z = k + 1
        x0, x1 = (0, 120) if z <= 10 else (20, 170) if z <= 20 else (40, 200)
        R = np.hypot(z, 3)
        fig, ax = plt.subplots(figsize=A4)
        channel_ax(ax, k, x0, x1)
        ax.set_title(f"Canal {z}  (z = {z} m)  -  onda S", fontsize=15)
        ax.axvline(T[k], color="green", lw=2.5)
        ax.text(T[k] + 1, ax.get_ylim()[1] * 0.85, f"t = {T[k]:g} ms", color="green", fontsize=14, weight="bold")
        fig.tight_layout(rect=(0, 0.1, 1, 1))
        fig.text(0.06, 0.04, f"t leído = {T[k]:g} ms        R = √({z}² + 3²) = {R:.2f} m        tz = {T[k]:g} · {z} / {R:.2f} = {TZ[k]:.1f} ms"
                 + ("   (no se usa: campo cercano)" if z <= 2 else "" if z <= 17 else "   (lectura relativa)"),
                 fontsize=13)
        fig.savefig(out / f"canal_{z:02d}.png", dpi=90); pdf.savefig(fig); plt.close(fig)
    # gráfica tz contra z
    L = ((3, 6), (6, 11), (11, 16), (17, 30))
    fig, ax = plt.subplots(figsize=(8, 9))
    ax.plot(TZ, Z, "o", color="#1f4e79", ms=6, label="tz leído")
    hs, vs_ = [], []
    for c1, c2 in L:
        m = (Z >= c1) & (Z <= c2); p = np.polyfit(Z[m], TZ[m], 1); v = 1000 / p[0]
        ax.plot(np.polyval(p, [c1, c2]), [c1, c2], color="#c0392b", lw=2)
        ax.text(np.polyval(p, (c1 + c2) / 2) + 4, (c1 + c2) / 2, f"Vs = {v:.0f} m/s", va="center", fontsize=11)
        vs_.append(v)
    h = np.diff([0, 6, 11, 16.5, 30]); vs30 = 30 / np.sum(h / np.array(vs_))
    ax.set_ylim(31, 0); ax.set_xlim(0, 200); ax.grid(alpha=.3)
    ax.set_xlabel("tz (ms)"); ax.set_ylabel("Profundidad z (m)")
    ax.set_title(f"Tiempo vertical contra profundidad  (Vs30 = {vs30:.0f} m/s)")
    fig.tight_layout(); fig.savefig(ROOT / "figuras" / "downhole_tz_vs_z.png", dpi=200); pdf.savefig(fig); plt.close(fig)
    print([round(v) for v in vs_], round(vs30))
print("ok")
