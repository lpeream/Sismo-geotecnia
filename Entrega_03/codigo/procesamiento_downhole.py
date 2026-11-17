"""Procesamiento del ensayo down hole (punto 1).

Lee los registros crudos de una copia de trabajo de recursos/DATOS DOWN HOLE
2026.xlsx (datos/downhole/, hojas VP y VS); solo lee, no escribe en el Excel,
obtiene los tiempos de llegada de las ondas P y S en cada canal y calcula los
perfiles de Vp y Vs, que se comparan con la gráfica de referencia del
enunciado.

Supuestos:
* 40 canales cada 1 m (canal k a k m de profundidad); se usan los 30 primeros.
* Fuente en superficie a 3 m del pozo (esquema de la hoja REFERENCIAS).
* Onda S: se usa (X - Y)/2, porque los golpes X y Y producen ondas S de
  polaridad opuesta y la resta las refuerza y cancela la onda P.
* Llegadas S (T_S_MS): canales 1 a 17 leídos a mano en las gráficas X-Y de
  cada canal, en el inicio de la onda S fuerte (X y Y en espejo con amplitud
  grande); las ondulaciones pequeñas previas son la onda P y no se leen.
  Canales 18 a 30: el inicio no es claro; se usan los retardos de la
  correlación cruzada entre canales vecinos, sumados a la lectura del canal 17.
* Los canales 1 y 2 no entran en el cálculo (campo cercano: a menos de una
  longitud de onda de la fuente y con la señal saturada).
* Llegadas P (T_P_MS): primer quiebre leído a mano en cada registro (ASTM
  D7400). El canal 1 de la hoja VP se lee de la caché del vínculo externo.
* Velocidad de cada tramo: inversa de la pendiente del tiempo vertical
  equivalente, t_v = t * z / sqrt(z^2 + d^2), contra la profundidad.
"""
from __future__ import annotations

import re
import zipfile
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.signal import hilbert

DELIVERY_ROOT = Path(__file__).resolve().parents[1]
# Se lee una copia de trabajo; el original de recursos/ no se usa ni se modifica.
RAW = DELIVERY_ROOT / "datos" / "downhole" / "DATOS DOWN HOLE 2026 - copia de trabajo.xlsx"
REFERENCE = DELIVERY_ROOT / "datos" / "ensayos_laboratorio" / "vs_downhole_leido_de_grafica.csv"
RESULTS = DELIVERY_ROOT / "datos" / "resultados"
FIGURES = DELIVERY_ROOT / "figuras"

N_CHANNELS = 40
DZ_M = 1.0
OFFSET_M = 3.0
MAX_DEPTH_M = 30.0
SEARCH_S = 0.25          # ventana de búsqueda desde el disparo
CC_WINDOW_S = 0.040      # ventana de correlación de la onda S directa
CC_PRE_S = 0.005
CC_MAX_LAG_S = 0.010
# Estratos del perfil: (techo, base, canales usados en la regresión). El cambio
# de pendiente de la onda S está entre los canales 16 y 17: el límite se toma en
# 16.5 m. En la arena las llegadas son claras y se usan tres rectas; en la
# arcilla (canales 17 a 30) la llegada no es clara y se usa una sola recta, para
# que los errores de lectura se compensen (con rectas cortas el error de la
# pendiente es de ±60 a ±90 m/s).
LAYERS = (
    (0.0, 6.0, (3, 6)), (6.0, 11.0, (6, 11)), (11.0, 16.5, (11, 16)),
    (16.5, 30.0, (17, 30)),
)
# Tiempos de llegada de la onda S (ms), canales 1 a 30 (ver figuras/lectura_downhole_S/).
# Canal 18: la correlación daba 129 ms (salto de un ciclo desde el canal 17); se
# ajustó a mano a 127 ms.
# Tiempos de llegada de la onda P (ms), canales 1 a 30, primer quiebre leído a mano.
T_P_MS = np.array([4.0, 8.0, 8.5, 9.0, 9.5, 10.0, 11.0, 11.5, 11.0, 11.0, 13.0, 13.0, 14.0,
                   14.0, 15.0, 15.5, 15.5, 16.0, 17.5, 17.5, 17.5, 19.0, 19.0, 20.0, 20.5,
                   20.5, 22.0, 22.5, 23.0, 24.0])
T_S_MS = np.array([18, 21, 36, 38.5, 41, 52.5, 54, 62.5, 63, 80, 80.5, 90.5, 99.5, 105, 104.5,
                   113.5, 119, 127, 129, 129, 133, 135.5, 136.5, 137, 144.5, 144.5, 147, 151.5,
                   153.5, 154.0])


def channel1_p_from_link_cache(n_samples):
    """Canal 1 de la onda P.

    En la hoja VP la columna B (canal 1) son fórmulas "=+'[1]VP (2)'!D9..." que
    apuntan a otro libro, por lo que no tienen valor guardado en la hoja. Excel
    conserva los valores de ese libro en xl/externalLinks/externalLink1.xml; se
    leen de ahí (fila 3 de VP = fila 9 de 'VP (2)').
    """
    with zipfile.ZipFile(RAW) as z:
        xml = z.read("xl/externalLinks/externalLink1.xml").decode("utf-8")
    names = re.findall(r'<sheetName val="([^"]+)"', xml)
    sheet_id = names.index("VP (2)")
    block = re.search(rf'<sheetData sheetId="{sheet_id}">(.*?)</sheetData>', xml, re.S).group(1)
    values = {int(r): float(v) for r, v in re.findall(r'<cell r="D(\d+)"><v>([^<]+)</v>', block)}
    return np.array([values.get(9 + i, np.nan) for i in range(n_samples)])


def read_records():
    vp = pd.read_excel(RAW, sheet_name="VP", header=None)
    vs = pd.read_excel(RAW, sheet_name="VS", header=None)
    t = vp.iloc[2:, 0].astype(float).to_numpy()
    p = np.array([vp.iloc[2:, 3 * k + 1].astype(float).to_numpy() for k in range(N_CHANNELS)])
    if np.isnan(p[0]).all():
        p[0] = channel1_p_from_link_cache(len(t))
    sx = np.array([vs.iloc[3:, 6 * k + 1].astype(float).to_numpy() for k in range(N_CHANNELS)])
    sy = np.array([vs.iloc[3:, 6 * k + 4].astype(float).to_numpy() for k in range(N_CHANNELS)])
    return t, p, sx, sy


def s_arrivals(t, s):
    """Retardos S por correlación cruzada; de aquí salen los canales 18 a 30 de T_S_MS."""
    dt = t[1] - t[0]
    n = int(SEARCH_S / dt)
    rough = []
    for x in s:
        env = np.abs(hilbert(x[:n] - np.median(x[:10])))
        rough.append(np.argmax(env > 0.3 * env.max()))
    rough = np.maximum.accumulate(np.array(rough))
    win, pre, max_lag = int(CC_WINDOW_S / dt), int(CC_PRE_S / dt), int(CC_MAX_LAG_S / dt)
    lags, scores = [], []
    for k in range(len(s) - 1):
        a0 = max(rough[k] - pre, 0)
        a = s[k][a0:a0 + win]
        best, best_lag = -2.0, 0
        for lag in range(max_lag + 1):
            b = s[k + 1][a0 + lag:a0 + lag + win]
            c = np.dot(a - a.mean(), b - b.mean()) / (np.std(a) * np.std(b) * len(a) + 1e-30)
            if c > best:
                best, best_lag = c, lag
        lags.append(best_lag * dt)
        scores.append(best)
    first = rough[0] * dt
    return np.r_[first, first + np.cumsum(lags)], np.r_[np.nan, scores]


def segment_velocities(z, t_vert):
    rows = []
    for top, bottom, (c1, c2) in LAYERS:
        m = (z >= c1) & (z <= c2)
        slope = np.polyfit(z[m], t_vert[m], 1)[0]
        rows.append((top, bottom, 1.0 / slope, c1, c2))
    return rows


def record_plot(t, traces, picks, title, filename, second=None, t_max=250.0):
    """Registros de todos los canales apilados por profundidad, con la llegada leída.

    Cada traza se centra en su valor inicial y se escala con su máximo dentro de
    la ventana mostrada; los registros vienen saturados (topes planos).
    """
    t_ms = t * 1000
    m = t_ms <= t_max
    fig, ax = plt.subplots(figsize=(11, 9), dpi=300)
    for k in range(N_CHANNELS):
        if np.isnan(traces[k]).all():
            continue
        sets = [(traces[k], "#1f4e79")] + ([(second[k], "#c0392b")] if second is not None else [])
        scale = max(np.nanmax(np.abs(x[m] - np.nanmedian(x[:5])), initial=0) for x, _ in sets) + 1e-30
        for x, color in sets:
            y = np.clip((x[m] - np.nanmedian(x[:5])) / scale, -1, 1)
            ax.plot(t_ms[m], (k + 1) - 0.45 * y, color=color, lw=0.5)
    for k in range(N_CHANNELS):
        if np.isnan(traces[k]).all():
            ax.text(0.02 * t_max, k + 1, f"Canal {k + 1}: sin datos en el Excel", va="center", fontsize=8,
                    color="gray")
    used = np.arange(1, len(picks) + 1)
    ax.plot(picks * 1000, used, "|", color="#2e8b57", ms=9, mew=2, label="Llegada leída")
    if second is not None:
        ax.plot([], [], color="#1f4e79", lw=1, label="Golpe X")
        ax.plot([], [], color="#c0392b", lw=1, label="Golpe Y")
    ax.set_xlim(0, t_max)
    ax.set_ylim(N_CHANNELS + 1, 0)
    ax.set_xlabel("Tiempo (ms)")
    ax.set_ylabel("Profundidad (m)")
    ax.set_title(title, fontsize=12, weight="bold")
    ax.grid(True, axis="x", alpha=0.3)
    ax.legend(loc="lower right", frameon=True, fontsize=9)
    fig.tight_layout()
    fig.savefig(FIGURES / filename, bbox_inches="tight")
    plt.close(fig)


def run_analysis() -> dict[str, pd.DataFrame]:
    t, p, sx, sy = read_records()
    z = np.arange(1, N_CHANNELS + 1) * DZ_M
    t_p = np.full(N_CHANNELS, np.nan)
    t_p[:len(T_P_MS)] = T_P_MS / 1000
    t_s = np.full(N_CHANNELS, np.nan)
    t_s[:len(T_S_MS)] = T_S_MS / 1000
    record_plot(t, p, t_p[:30], "Down hole - onda P", "downhole_registros_p.png", t_max=1000.0)
    record_plot(t, sx, T_S_MS / 1000, "Down hole - onda S", "downhole_registros_s.png", second=sy,
                t_max=1000.0)
    factor = z / np.hypot(z, OFFSET_M)
    arrivals = pd.DataFrame({
        "canal": np.arange(1, N_CHANNELS + 1), "z_m": z,
        "t_P_ms": t_p * 1000, "t_S_ms": t_s * 1000,
        "tv_P_ms": t_p * factor * 1000, "tv_S_ms": t_s * factor * 1000,
        "lectura_S": ["visual"] * 17 + ["relativa al canal 17"] * 13 + [""] * (N_CHANNELS - 30),
    })
    use = arrivals[arrivals["z_m"] <= MAX_DEPTH_M]
    zz = use["z_m"].to_numpy()
    vs_rows = segment_velocities(zz, use["tv_S_ms"].to_numpy() / 1000)
    vp_rows = segment_velocities(zz, use["tv_P_ms"].to_numpy() / 1000)
    profile = pd.DataFrame([
        {"desde_m": a, "hasta_m": b, "canales": f"{c1}-{c2}", "Vs_m_s": vs, "Vp_m_s": vp,
         "nu": (0.5 * (vp / vs) ** 2 - 1) / ((vp / vs) ** 2 - 1)}
        for (a, b, vs, c1, c2), (_, _, vp, _, _) in zip(vs_rows, vp_rows)
    ])
    travel = ((profile["hasta_m"] - profile["desde_m"]) / profile["Vs_m_s"]).sum()
    vs30 = 30.0 / travel

    ref = pd.read_csv(REFERENCE, encoding="utf-8-sig")
    RESULTS.mkdir(parents=True, exist_ok=True)
    arrivals.to_csv(RESULTS / "downhole_tiempos_llegada.csv", index=False, encoding="utf-8-sig")
    profile.to_csv(RESULTS / "downhole_perfil_velocidades.csv", index=False, encoding="utf-8-sig")

    # Figura: tiempo vertical contra profundidad y perfil de Vs contra la referencia.
    plt.rcParams.update({"font.size": 11})
    fig, (ax_t, ax_v) = plt.subplots(1, 2, figsize=(11, 7), dpi=300, sharey=True)
    excl = zz < LAYERS[0][2][0]
    ax_t.plot(use["tv_S_ms"][~excl], zz[~excl], "o", color="#1f4e79", ms=4, label="Onda S")
    ax_t.plot(use["tv_S_ms"][excl], zz[excl], "o", mfc="white", color="#1f4e79", ms=4,
              label="Onda S (no usada, campo cercano)")
    ax_t.plot(use["tv_P_ms"], zz, "s", color="#c08a2b", ms=4, label="Onda P")
    for (_, _, _, c1, c2) in vs_rows:
        m = (zz >= c1) & (zz <= c2)
        fit = np.polyfit(zz[m], use["tv_S_ms"].to_numpy()[m], 1)
        ax_t.plot(np.polyval(fit, [c1, c2]), [c1, c2], color="#1f4e79", lw=1.5)
    ax_t.set_xlabel("Tiempo vertical (ms)")
    ax_t.set_ylabel("Profundidad (m)")
    ax_t.set_title("Tiempos de llegada", fontsize=11)
    ax_t.legend(frameon=False)
    ax_t.grid(True, alpha=0.3)

    def steps(df, col, top="desde_m", bottom="hasta_m"):
        d = [v for _, r in df.iterrows() for v in (r[top], r[bottom])]
        s = [v for _, r in df.iterrows() for v in (r[col], r[col])]
        return s, d

    s_ref, d_ref = steps(ref, "Vs_m_s")
    s_own, d_own = steps(profile, "Vs_m_s")
    ax_v.plot(s_ref, d_ref, color="gray", lw=2, ls="--", label="Referencia del enunciado")
    ax_v.plot(s_own, d_own, color="#1f4e79", lw=2.5, label=f"Este trabajo (Vs30 ≈ {vs30:.0f} m/s)")
    ax_v.set_xlim(0, 500)
    ax_v.set_xlabel("Vs (m/s)")
    ax_v.set_title("Perfil de Vs", fontsize=11)
    ax_v.legend(frameon=False, loc="lower left")
    ax_v.grid(True, alpha=0.3)
    ax_v.set_ylim(MAX_DEPTH_M, 0)
    fig.tight_layout()
    fig.savefig(FIGURES / "downhole_tiempos_y_perfil_vs.png", bbox_inches="tight")
    plt.close(fig)

    print(arrivals[arrivals.z_m <= MAX_DEPTH_M].round(1).to_string(index=False))
    print()
    print(profile.round(2).to_string(index=False))
    print(f"\nVs30 propio = {vs30:.1f} m/s")
    return {"llegadas": arrivals, "perfil": profile}


if __name__ == "__main__":
    run_analysis()
