"""Escribe en el informe las secciones de los puntos 1 y 2 (texto y tablas).

Lee los ensayos (datos/ensayos_laboratorio/) y los resultados de los scripts
(datos/resultados/) y reemplaza en informe/entrega_03.tex todo lo que va desde
"\\subsection{Ensayos disponibles}" hasta "\\section{Punto 3". Correr después de
procesamiento_downhole.py, perfil_estratigrafico.py y curvas_degradacion.py.
"""
from pathlib import Path
import pandas as pd

R = Path(__file__).resolve().parents[1]   # carpeta Entrega_03
LAB = R / "datos" / "ensayos_laboratorio"
RES = R / "datos" / "resultados"
TEX = R / "informe" / "entrega_03.tex"
NL = " \\\\"          # fin de fila LaTeX


def f(x, d=1):
    s = f"{x:,.{d}f}".replace(",", "\\,")
    return s


n_tab = [0]


def table(title, label, cols, header, rows, size="\\footnotesize"):
    """Tabla con cuadrícula completa y encabezados en negrita (estilo de la referencia)."""
    n_tab[0] += 1
    spec = "|" + "|".join("c" if c in "rc" else (r">{\raggedright\arraybackslash}" + c if c.startswith("p") else c) for c in _split_cols(cols)) + "|"
    head = " & ".join("\\thead{" + h.strip() + "}" for h in header.split(" & "))
    out = ["\\begin{table}[H]", "\\centering",
           f"\\customtablecaption{{Tabla {n_tab[0]}}}{{{title}}}", f"\\label{{{label}}}", size,
           "\\renewcommand{\\arraystretch}{1.25}",
           f"\\begin{{tabular}}{{{spec}}}", "\\hline", head + NL, "\\hline"]
    out += [r + NL + " \\hline" for r in rows]
    out += ["\\end{tabular}", "\\end{table}", ""]
    return "\n".join(out)


def _split_cols(cols):
    """Separa una especificación como 'lrp{2cm}c' en columnas."""
    out, i = [], 0
    while i < len(cols):
        if cols[i] == "p":
            j = cols.index("}", i)
            out.append(cols[i:j + 1]); i = j + 1
        else:
            out.append(cols[i]); i += 1
    return out


# ---------------- datos ----------------
cl = pd.read_csv(LAB / "clasificacion_muestras.csv", encoding="utf-8-sig")
rc = pd.read_csv(LAB / "columna_resonante.csv", encoding="utf-8-sig")
tx = pd.read_csv(LAB / "triaxial_ciclico_5m.csv", encoding="utf-8-sig")
cu = pd.read_csv(LAB / "triaxial_cu_parametros.csv", encoding="utf-8-sig")
vel = pd.read_csv(RES / "downhole_perfil_velocidades.csv", encoding="utf-8-sig")
tram = pd.read_csv(RES / "perfil_tramos_vs_gmax.csv", encoding="utf-8-sig")
par = pd.read_csv(RES / "parametros_curvas_degradacion.csv", encoding="utf-8-sig")


def cons_il(il):
    return "Semisólido o sólido" if il < 0 else "Firme" if il < 0.5 else "Blando" if il < 1 else "Líquido"


def cons_ic(ic):
    return ("Duro" if ic > 1 else "Rígido" if ic >= 0.75 else "Medio" if ic >= 0.5
            else "Blando" if ic >= 0.25 else "Muy blando")


def plast(ip):
    return "Poco plástico" if ip < 7 else "Medianamente plástico" if ip <= 17 else "Altamente plástico"


blocks = []

# ---------------- 1.1 Ensayos ----------------
blocks.append(r"""\subsection{Ensayos disponibles}

Se cuenta con un sondeo (P-1) de 30 m con los ensayos de la
Tabla~\ref{tab:ensayos}.

""")
blocks.append(table("Ensayos de campo y laboratorio del sondeo P-1.", "tab:ensayos", "lp{2.6cm}p{4.2cm}",
                    "Ensayo & Profundidades\\\\(m) & Uso",
                    ["Down hole ($V_p$ y $V_s$, 40 canales) & 1 a 40 & $V_s$, $V_p$, $\\nu$ y $G_{\\max}$",
                     "Humedad, granulometría y límites & 3, 6, 12, 16, 22, 26 y 30 & Clasificación USCS",
                     "Triaxial estático CU & 4 y 28 & $\\phi'$, $c'$ y $\\gamma_h$",
                     "Triaxial cíclico (ASTM D3999) & 5 & Curva de la arena",
                     "Columna resonante (ASTM D4015) & 10 y 20 & Curvas de arena y arcilla"]))

# ---------------- 1.2 Clasificación ----------------
rows = [f"{r.muestra} & {r.profundidad_m:.0f} & {r.humedad_pct:.1f} & {r.grava_pct:.1f} & {r.arena_pct:.1f} & "
        f"{r.finos_pct:.1f} & {r.LL_pct:.0f} & {r.LP_pct:.0f} & {r.IP_pct:.0f} & {r.USCS_informado}"
        for r in cl.itertuples()]
blocks.append(r"""\subsection{Clasificación USCS}

Se ven dos estratos claramente diferentes (Tabla~\ref{tab:clasificacion}):
arriba una arena con grava y limo (SP-SM) y abajo una arcilla. La arcilla es de
baja plasticidad (CL) y se vuelve más plástica a medida que se baja (LL de 26
a 50\%); la muestra de 30 m ya clasifica como CH.

""")
blocks.append(table("Clasificación de las muestras del sondeo P-1.", "tab:clasificacion", "lrrrrrrrrl",
                    "Muestra & Prof.\\\\(m) & $w$\\\\(\\%) & Grava\\\\(\\%) & Arena\\\\(\\%) & Finos\\\\(\\%) & LL\\\\(\\%) & LP\\\\(\\%) & IP\\\\(\\%) & USCS",
                    rows))
blocks.append(r"""\pendiente{Revisar con el profesor: la muestra P1-M12 tiene la misma
granulometría y los mismos límites que P1-M3, y en ambas la grava (47.1\%)
supera a la arena (43.3\%), por lo que la USCS daría GP-GM o GW-GM.}

""")

# ---------------- 1.3 Atterberg ----------------
rows = []
for r in cl.itertuples():
    il = (r.humedad_pct - r.LP_pct) / r.IP_pct
    ic = (r.LL_pct - r.humedad_pct) / r.IP_pct
    rows.append(f"{r.muestra} & {r.profundidad_m:.0f} & {r.humedad_pct:.1f} & {r.IP_pct:.0f} & "
                f"{r.USCS_informado} & {il:.2f} & {ic:.2f} & {cons_il(il)} & {cons_ic(ic)} & {plast(r.IP_pct)}")
blocks.append(r"""\subsection{Consistencia a partir de los límites de Atterberg}

Para no quedarnos solo con el IP, se calculan también el índice de liquidez,
$I_L=(w-LP)/IP$, y el de consistencia, $I_C=(LL-w)/IP$
(Tabla~\ref{tab:atterberg}).

""")
blocks.append(table("Índices de liquidez y consistencia.", "tab:atterberg", "lrrrlrrp{2.1cm}p{1.5cm}p{2.3cm}",
                    "Muestra & Prof.\\\\(m) & $w$\\\\(\\%) & IP\\\\(\\%) & USCS & $I_L$ & $I_C$ & Consist.\\\\según $I_L$ & Consist.\\\\según $I_C$ & Plasticidad",
                    rows))
blocks.append(r"""\begin{itemize}
    \item \textbf{Estrato 1 (arena):} es poco plástica y de firme a dura. Como
    es un suelo granular, lo que manda en su comportamiento es la fricción y
    no los finos.
    \item \textbf{Estrato 2 (arcilla):} la plasticidad sube con la
    profundidad, de media a alta. Entre 16 y 22 m es blanda; a 22 m la humedad
    (36.6\%) incluso supera el límite líquido. Hacia el fondo se vuelve más
    firme.
\end{itemize}

""")

# ---------------- 1.4 Triaxial CU ----------------
rows = []
for r, e in zip(cu.itertuples(), ("E1 (arena)", "E2 (arcilla)")):
    rows.append(f"{e} & {r.profundidad_m:.0f} & {r.peso_unitario_kN_m3:.1f} & {r.phi_efectivo_grados:.0f} & {r.c_efectiva_kPa:.0f}")
blocks.append(r"""\subsection{Resistencia a partir del triaxial estático CU}

Se hicieron dos triaxiales CU, uno en cada estrato (Tabla~\ref{tab:cu}). El
criterio de falla fue el máximo esfuerzo desviador.

""")
blocks.append(table("Parámetros de resistencia del triaxial CU.", "tab:cu", "lrrrr",
                    "Estrato & Prof.\\\\(m) & $\\gamma_h$\\\\(kN/m$^3$) & $\\phi'$\\\\(°) & $c'$\\\\(kPa)", rows))
blocks.append(r"""La arena no tiene cohesión, como es de esperar en un suelo granular. La
arcilla tiene cohesión y casi nada de fricción, así que su resistencia es
prácticamente la no drenada.

""")

# ---------------- 1.5 Parámetros dinámicos ----------------
blocks.append(r"""\subsection{Parámetros dinámicos}

\subsubsection{Ensayo Down Hole}

El down hole tiene 40 geófonos, uno por metro, y la fuente está en superficie
a 3 m del pozo. Para la onda P se golpea en vertical y para la onda S se dan
dos golpes horizontales en sentidos contrarios (X y Y). Las Figuras~\ref{fig:dh_p}
y~\ref{fig:dh_s} muestran los registros, con la llegada leída en cada canal.

\begin{figure}[H]
    \centering
    \includegraphics[width=0.9\linewidth]{../figuras/downhole_registros_p.png}
    \caption{Down hole: onda P.}
    \label{fig:dh_p}
\end{figure}

\begin{figure}[H]
    \centering
    \includegraphics[width=0.9\linewidth]{../figuras/downhole_registros_s.png}
    \caption{Down hole: onda S.}
    \label{fig:dh_s}
\end{figure}

La llegada de la onda P se leyó en el primer quiebre de cada registro
\parencite{astmD7400}, y la de la onda S donde los golpes X y Y se separan en
espejo. Como la fuente no está sobre el pozo, cada tiempo se pasa a vertical
con $t_z=t\,z/\sqrt{z^2+3^2}$, y la velocidad sale de la pendiente de $t_z$
contra la profundidad. En la arena las llegadas son claras y se usan tres
rectas; en la arcilla la llegada no se ve tan clara canal por canal, así que
se usa una sola recta con los 14 canales para que los errores de lectura se
compensen. Los dos primeros metros no se usan porque están muy cerca de la
fuente. Con las velocidades se calculan el coeficiente de Poisson
$\nu$ y el módulo cortante $G_{\max}=\rho V_s^2$ de cada estrato
(Tabla~\ref{tab:perfil}). Para $\rho$ se usa el peso unitario promedio de las
muestras de cada estrato (Tabla~\ref{tab:gamma}). En los triaxiales CU el
peso unitario viene directo en los resultados; en el triaxial cíclico y en la
columna resonante el ensayo da la densidad húmeda, y se pasa a peso unitario
con $\gamma_h=\rho_h\,g$. Por ejemplo, a 5 m: $\gamma_h=1.891\ \text{g/cm}^3
\times9.81=18.55$ kN/m$^3$.

""")
# Peso unitario de las muestras y promedio por estrato.
din = pd.read_csv(LAB / "propiedades_muestras_dinamicas.csv", encoding="utf-8-sig")
g_rows = [("E1", 4.0, "Triaxial CU", float(cu.peso_unitario_kN_m3.iloc[0]))]
g_rows += [("E1" if r.profundidad_m < 16.5 else "E2", r.profundidad_m, r.ensayo.split(" ASTM")[0],
            r.densidad_humeda_g_cm3 * 9.81) for r in din.itertuples()]
g_rows += [("E2", 28.0, "Triaxial CU", float(cu.peso_unitario_kN_m3.iloc[1]))]
g_rows.sort(key=lambda x: x[1])
rows = []
for est in ("E1", "E2"):
    sub = [g for g in g_rows if g[0] == est]
    for e, z, ens, g in sub:
        rows.append(f"{e[1]} & {z:g} & {ens} & {g:.2f}")
    rows.append(f"\\textbf{{Promedio estrato {est[1]}}} & & & \\textbf{{{sum(x[3] for x in sub) / len(sub):.1f}}}")
blocks.append(table("Peso unitario de las muestras.", "tab:gamma", "cclc",
                    "Estrato & Prof.\\\\(m) & Ensayo & $\\gamma_h$\\\\(kN/m$^3$)",
                    rows))
# Valores por estrato: promedio por tiempo de viaje de las capas del down hole.
rows = []
for est, sub_t in tram.groupby("estrato"):
    h = sub_t.hasta_m - sub_t.desde_m
    sub_v = vel[(vel.desde_m >= sub_t.desde_m.min()) & (vel.hasta_m <= sub_t.hasta_m.max())]
    vs_e = h.sum() / (h / sub_t.Vs_m_s).sum()
    vp_e = h.sum() / ((sub_v.hasta_m - sub_v.desde_m) / sub_v.Vp_m_s).sum()
    r2 = (vp_e / vs_e) ** 2
    nu_e = (r2 - 2) / (2 * (r2 - 1))
    gam = sub_t.peso_unitario_kN_m3.iloc[0]
    rho = gam / 9.81 * 1000
    rows.append(f"{est[1]} & {h.sum():g} & {f(vp_e, 0)} & {vs_e:.0f} & {nu_e:.3f} & {gam:.1f} & {f(rho, 0)} & {rho * vs_e**2 / 1e6:.1f}")
    if est == "E1":
        vs_e1 = vs_e
    else:
        vs_e2 = vs_e
blocks.append(table("Propiedades de los estratos.", "tab:perfil", "cccccccc",
                    "Estrato & Espesor\\\\(m) & $V_p$\\\\(m/s) & $V_s$\\\\(m/s) & $\\nu$ & $\\gamma_h$\\\\(kN/m$^3$) & $\\rho$\\\\(kg/m$^3$) & $G_{\\max}$\\\\(MPa)",
                    rows))
blocks.append(r"""Como $\nu$ está muy cerca de 0.5 y $V_p$ es alta, el suelo está saturado y
el nivel freático está cerca de la superficie. El perfil de $V_s$ obtenido es muy
parecido al de la gráfica de referencia del enunciado
(Figura~\ref{fig:perfil}).

""")

# Triaxial cíclico
rows = []
for i, r in enumerate(tx.itertuples(), 1):
    tau = r.G_kPa * r.gamma_pct / 100
    rows.append(f"{i} & {r.gamma_pct:g} & {r.D_pct:.2f} & {f(r.E_kPa, 1)} & {f(r.G_kPa, 1)} & {tau:.1f}")
blocks.append(r"""\subsubsection{Triaxial cíclico}

Se hizo a una muestra de 5 m de profundidad (P1-M5), que está en el estrato 1,
con un confinamiento de 200 kPa (Tabla~\ref{tab:tc}).

""")
blocks.append(table("Resultados del triaxial cíclico, muestra P1-M5 (5 m).", "tab:tc", "rrrrrr",
                    "Ensayo & $\\gamma$\\\\(\\%) & $D$\\\\(\\%) & $E$\\\\(kPa) & $G$\\\\(kPa) & $\\tau$\\\\(kPa)", rows))

# Columna resonante
blocks.append(r"""\subsubsection{Columna resonante}

Se hicieron dos ensayos, uno en cada estrato: a 10 m (P1-M10) con 200 kPa
de confinamiento y a 20 m (P1-M20) con 400 kPa (Tablas~\ref{tab:cr10}
y~\ref{tab:cr20}).

""")
for z, lab, est in ((10.0, "tab:cr10", "P1-M10 (10 m)"), (20.0, "tab:cr20", "P1-M20 (20 m)")):
    sub = rc[rc.profundidad_m == z]
    rows = [f"{i} & {r.frecuencia_Hz:.0f} & {r.torque_Nm:.2f} & {r.gamma_pct:g} & {r.G_MPa:.2f} & {r.D_pct:.2f}"
            for i, r in enumerate(sub.itertuples(), 1)]
    blocks.append(table(f"Resultados de la columna resonante, muestra {est}.", lab, "rrrrrr",
                        "Ensayo & $f_r$\\\\(Hz) & Torque\\\\(N$\\cdot$m) & $\\gamma$\\\\(\\%) & $G$\\\\(MPa) & $D$\\\\(\\%)", rows))

# ---------------- 1.6 Perfil adoptado ----------------
blocks.append(r"""\subsection{Perfil adoptado}

El perfil tiene dos estratos (Figura~\ref{fig:perfil}): una arena con grava y
limo de 16.5 m y debajo una arcilla más rígida de 13.5 m. La muestra de 30 m
sale CH, pero está justo en el límite con CL (LL $=50\%$) y su $V_s$ es igual
al del resto de la arcilla, así que se deja en el mismo estrato.

La velocidad media de onda cortante en los 30 m superiores es
\begin{equation*}
\bar{V}_s=\frac{\sum d_i}{\sum d_i/V_{si}}
=\frac{30}{\frac{16.5}{""" + f"{vs_e1:.1f}" + r"""}+\frac{13.5}{""" + f"{vs_e2:.1f}" + r"""}}
\approx207\ \text{m/s},
\end{equation*}
por lo que el perfil es tipo D según la NSR-10 (entre 180 y 360 m/s).

\begin{figure}[H]
    \centering
    \includegraphics[width=0.72\linewidth]{../figuras/perfil_estratigrafico_p1.png}
    \caption{Perfil estratigráfico.}
    \label{fig:perfil}
\end{figure}

""")

# ---------------- Punto 2 ----------------
blocks.append(r"""\section{Punto 2: curvas de degradación}

Las curvas se arman con los datos del down hole, la columna resonante y el
triaxial cíclico. Del down hole sale $G_{\max}=\rho V_s^2$, que sirve para
normalizar los $G$ de los otros ensayos. Se usa el $G_{\max}$ de cada
estrato de la Tabla~\ref{tab:perfil} (Tablas~\ref{tab:dinE1}
y~\ref{tab:dinE2}).

""")
def gmax_estrato(est):
    sub = tram[tram.estrato == est]
    h = sub.hasta_m - sub.desde_m
    vs_eq = h.sum() / (h / sub.Vs_m_s).sum()
    return sub.peso_unitario_kN_m3.iloc[0] / 9.81 * vs_eq**2 / 1000


dh1, dh2 = gmax_estrato("E1"), gmax_estrato("E2")
rows = [f"Down hole (estrato 1) & 1 & $\\approx0$ & -- & {dh1:.2f} & 1.000"]
for i, r in enumerate(rc[rc.profundidad_m == 10].itertuples(), 1):
    rows.append(f"Columna resonante 10 m & {i} & {r.gamma_pct:g} & {r.D_pct:.2f} & {r.G_MPa:.2f} & {r.G_MPa / dh1:.3f}")
for i, r in enumerate(tx.itertuples(), 1):
    rows.append(f"Triaxial cíclico 5 m & {i} & {r.gamma_pct:g} & {r.D_pct:.2f} & {r.G_kPa / 1000:.2f} & {r.G_kPa / 1000 / dh1:.3f}")
blocks.append(table(f"Parámetros dinámicos del estrato E1 (arena), $G_{{\\max}}={dh1:.2f}$ MPa.", "tab:dinE1", "lrrrrr",
                    "Ensayo & No. & $\\gamma$\\\\(\\%) & $D$\\\\(\\%) & $G$\\\\(MPa) & $G/G_{\\max}$", rows))
sub = rc[rc.profundidad_m == 20]
g0 = sub.G_MPa.iloc[0]
rows = [f"Down hole (estrato 2) & 1 & $\\approx0$ & -- & {dh2:.2f} & 1.000"]
for i, r in enumerate(sub.itertuples(), 1):
    rows.append(f"Columna resonante 20 m & {i} & {r.gamma_pct:g} & {r.D_pct:.2f} & {r.G_MPa:.2f} & {r.G_MPa / dh2:.3f}")
blocks.append(table(f"Parámetros dinámicos del estrato E2 (arcilla), $G_{{\\max}}={dh2:.2f}$ MPa.", "tab:dinE2", "lrrrrr",
                    "Ensayo & No. & $\\gamma$\\\\(\\%) & $D$\\\\(\\%) & $G$\\\\(MPa) & $G/G_{\\max}$", rows))
blocks.append(r"""En los dos estratos el primer punto de la columna resonante queda por debajo
del $G_{\max}$ del down hole: 0.93 en la arena (40.1 contra 43.1 MPa) y 0.96
en la arcilla (264.4 contra 275.7 MPa). Es normal que la muestra salga algo
más blanda que el suelo en el sitio, porque se altera al sacarla. Por eso la
forma de la curva se toma de los ensayos (cada $G$ dividido por el primero de
la columna resonante) y el $G_{\max}$ se toma del down hole. En la
Figura~\ref{fig:curvas} los puntos de laboratorio se dibujan con esa forma, por
lo que caen sobre la curva.

""")

e1, e2 = par.iloc[0], par.iloc[1]
blocks.append(r"""Para la arena hay 13 puntos que van desde deformaciones muy pequeñas hasta
grandes, así que la curva se ajusta directamente a los ensayos. Para la
arcilla solo hay columna resonante, que llega hasta $\gamma=0.05\%$, así que
hay que completar la curva con un modelo de la literatura. En los dos casos se
usa el modelo hiperbólico
\begin{equation*}
\frac{G}{G_{\max}}=\frac{1}{1+\left(\gamma/\gamma_r\right)^{a}},
\end{equation*}
donde $\gamma_r$ es la deformación de referencia (en la que $G/G_{\max}=0.5$)
y $a$ la curvatura:
\begin{itemize}
    \item \textbf{E1 (arena):} $\gamma_r$ y $a$ se ajustan a la columna
    resonante de 10 m y al triaxial cíclico de 5 m. El amortiguamiento se
    ajusta con $D=D_{\min}+(D_{\max}-D_{\min})(1-G/G_{\max})^{n}$, con
    $D_{\max}=17.0\%$ y $n=0.98$.
    \item \textbf{E2 (arcilla):} se ajusta $\gamma_r$ a la columna resonante
    de 20 m con la curvatura del modelo de \textcite{darendeli2001}
    ($a=0.919$), que extrapola la curva a deformaciones grandes; el
    amortiguamiento sigue la expresión de ese mismo modelo.
\end{itemize}

""")
rows = [f"E1 & Columna resonante 10 m y triaxial cíclico 5 m & {e1.IP:.1f} & {e1.gamma_r_pct:.2f} & {e1.a:.2f} & {e1.D_min_pct:.2f}",
        f"E2 & Columna resonante 20 m (forma de Darendeli) & {e2.IP:.1f} & {e2.gamma_r_pct:.2f} & {e2.a:.3f} & {e2.D_min_pct:.2f}"]
blocks.append(table("Parámetros de las curvas de degradación.", "tab:curvas", "lp{5.5cm}rrrr",
                    "Estrato & Origen & IP\\\\(\\%) & $\\gamma_r$\\\\(\\%) & $a$ & $D_{\\min}$\\\\(\\%)", rows))
blocks.append(r"""\begin{figure}[H]
    \centering
    \includegraphics[width=\linewidth]{../figuras/curvas_degradacion_estratos.png}
    \caption{Curvas de degradación y amortiguamiento por estrato.}
    \label{fig:curvas}
\end{figure}

La curva de la arena se ajusta bien a los puntos de los ensayos
($R^2=0.9999$). Los dos suelos se degradan menos de lo que predice el modelo
de Darendeli; esto puede deberse a que las muestras se ensayaron con más
confinamiento que el que tienen en campo, por lo que en el sitio las curvas
podrían caer un poco antes.
\pendiente{Nota para el punto 5 (DEEPSOIL): la curva de la arcilla (E2)
solo tiene datos hasta $\gamma=0.05\%$ (columna resonante de 20 m); el resto
(línea punteada) se completó por ahora con la forma del modelo de Darendeli
($a=0.919$, $\gamma_r=0.65\%$) ajustada a esos puntos. La idea es reemplazarla,
como en el trabajo de referencia, por la curva de Vucetic y Dobry (1991) que
trae DEEPSOIL: probar IP $=15$ y IP $=30$ (la arcilla tiene IP $\approx18$),
compararlas con los 5 puntos de la columna resonante (Tabla~\ref{tab:dinE2}) y
escoger la que mejor se ajuste. Al cambiarla hay que actualizar la
Figura~\ref{fig:curvas}, la Tabla~\ref{tab:curvas} y este texto. Revisar
también $D_{\min}$ de la arcilla: el ajuste dio 0\% porque el ensayo mide
0.02\% a $\gamma=0.002\%$, muy bajo para DEEPSOIL (valores usuales de 0.5 a
2\%; la referencia usó 5\% en todas las capas). La arena (E1) no necesita
cambios: su curva sale completa de los ensayos.}

""")

s = TEX.read_text(encoding="utf-8")
a = s.index("\\subsection{Ensayos disponibles}")
b = s.index("\\section{Punto 3")
s = s[:a] + "".join(blocks) + s[b:]
TEX.write_text(s, encoding="utf-8")
print("tablas:", n_tab[0])
