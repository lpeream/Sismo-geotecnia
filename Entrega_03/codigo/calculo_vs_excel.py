"""Excel con el cálculo de Vp y Vs del down hole, con fórmulas editables.

Genera datos/downhole/calculo_vs_downhole.xlsx a partir de los tiempos de
llegada (datos/resultados/downhole_tiempos_llegada.csv). Las celdas amarillas
(tiempos, canales de cada capa, peso unitario) se pueden cambiar en Excel.
"""
from pathlib import Path

import openpyxl, pandas as pd
from openpyxl.styles import Font, PatternFill, Alignment, Border, Side
from openpyxl.chart import ScatterChart, Reference, Series
ROOT = Path(__file__).resolve().parents[1]   # carpeta Entrega_03
F = ROOT / "datos" / "downhole" / "calculo_vs_downhole.xlsx"
arr = pd.read_csv(ROOT / "datos" / "resultados" / "downhole_tiempos_llegada.csv", encoding="utf-8-sig")
arr = arr[arr.z_m <= 30]
wb = openpyxl.Workbook(); ws = wb.active; ws.title = "CALCULO_VS"
ref = pd.read_csv(ROOT / "datos" / "ensayos_laboratorio" / "vs_downhole_leido_de_grafica.csv", encoding="utf-8-sig")
b = Font(bold=True); yel = PatternFill("solid", fgColor="FFF2CC"); blu = PatternFill("solid", fgColor="DDEBF7")
thin = Side(style="thin", color="999999"); box = Border(thin, thin, thin, thin)
ws["A1"] = "Cálculo de Vp y Vs del down hole"; ws["A1"].font = Font(bold=True, size=13)
ws["A2"] = "Amarillo = datos que se pueden cambiar (tiempos leídos). Lo demás son fórmulas."
ws["A3"] = "x fuente-pozo (m)"; ws["C3"] = 3; ws["C3"].fill = yel
ws["E3"] = "g (m/s²)"; ws["F3"] = 9.81
H = ["Canal", "z (m)", "t P (ms)", "t S (ms)", "R = √(z²+x²) (m)", "tv P = t·z/R (ms)",
     "tv S = t·z/R (ms)", "Vp prom = R/t (m/s)", "Vs prom = R/t (m/s)", "Vp intervalo = Δz/Δtv (m/s)", "Vs intervalo = Δz/Δtv (m/s)", "Vs referencia (m/s)", "Diferencia Vs (%)"]
for j, h in enumerate(H, 1):
    c = ws.cell(5, j, h); c.font = b; c.fill = blu; c.border = box
    c.alignment = Alignment(wrap_text=True, horizontal="center")
for i, r in enumerate(arr.itertuples(), 6):
    ws.cell(i, 1, int(r.canal)); ws.cell(i, 2, float(r.z_m))
    ws.cell(i, 3, round(float(r.t_P_ms), 2)).fill = yel
    ws.cell(i, 4, round(float(r.t_S_ms), 2)).fill = yel
    ws.cell(i, 5, f"=SQRT(B{i}^2+$C$3^2)")
    ws.cell(i, 6, f"=C{i}*B{i}/E{i}"); ws.cell(i, 7, f"=D{i}*B{i}/E{i}")
    ws.cell(i, 8, f"=E{i}/(C{i}/1000)"); ws.cell(i, 9, f"=E{i}/(D{i}/1000)")
    if i == 6:
        ws.cell(i, 10, f"=B{i}/(F{i}/1000)"); ws.cell(i, 11, f"=B{i}/(G{i}/1000)")
    else:
        ws.cell(i, 10, f'=IFERROR((B{i}-B{i-1})/((F{i}-F{i-1})/1000),"Δt=0")')
        ws.cell(i, 11, f'=IFERROR((B{i}-B{i-1})/((G{i}-G{i-1})/1000),"Δt=0")')
    zr = float(r.z_m) - 0.5
    ws.cell(i, 12, float(ref[(ref.desde_m <= zr) & (zr < ref.hasta_m)].Vs_m_s.iloc[0])).font = b
    ws.cell(i, 13, f'=IFERROR((K{i}-L{i})/L{i}*100,"")').number_format = "0"
    for j in range(1, 14):
        ws.cell(i, j).border = box
        if j >= 5: ws.cell(i, j).number_format = "0.00" if j <= 7 else "0"
        if j == 11: ws.cell(i, j).font = b
last = 5 + len(arr)
r0 = last + 3
ws.cell(r0 - 1, 1, "Velocidad por capa: V = 1 / pendiente de tv contra z (canal inicial y final editables)").font = b
H2 = ["Capa", "Desde (m)", "Hasta (m)", "Canal ini", "Canal fin", "Vs (m/s)", "Vp (m/s)", "ν",
      "h/Vs (s)", "γ (kN/m³)", "Gmax (MPa)", "Vs referencia (m/s)"]
for j, h in enumerate(H2, 1):
    c = ws.cell(r0, j, h); c.font = b; c.fill = blu; c.border = box
    c.alignment = Alignment(wrap_text=True, horizontal="center")
L = [(0, 6, 3, 6, 18.6, "143"), (6, 11, 6, 11, 18.6, "145–150"), (11, 16.5, 11, 16, 18.6, "150–160"), (16.5, 30, 17, 30, 18.7, "362–390")]
G, Bz, C = f"$G$6:$G${last}", f"$B$6:$B${last}", f"$F$6:$F${last}"
for k, (a, bb, c1, c2, gam, ref) in enumerate(L, 1):
    i = r0 + k
    ws.cell(i, 1, k); ws.cell(i, 2, a); ws.cell(i, 3, bb)
    ws.cell(i, 4, c1).fill = yel; ws.cell(i, 5, c2).fill = yel
    rg = lambda R: f"INDEX({R},D{i}):INDEX({R},E{i})"
    ws.cell(i, 6, f"=1000/SLOPE({rg(G)},{rg(Bz)})")
    ws.cell(i, 7, f"=1000/SLOPE({rg(C)},{rg(Bz)})")
    ws.cell(i, 8, f"=((G{i}/F{i})^2-2)/(2*((G{i}/F{i})^2-1))")
    ws.cell(i, 9, f"=(C{i}-B{i})/F{i}")
    ws.cell(i, 10, gam).fill = yel
    ws.cell(i, 11, f"=J{i}/$F$3*F{i}^2/1000")
    ws.cell(i, 12, ref)
    for j, fmt in zip(range(1, 13), ["0", "0.0", "0.0", "0", "0", "0", "0", "0.00", "0.0000", "0.0", "0", "@"]):
        ws.cell(i, j).number_format = fmt; ws.cell(i, j).border = box
rv = r0 + len(L) + 2
ws.cell(rv, 1, "Vs30 = 30 / Σ(h/Vs)").font = b
ws.cell(rv, 6, f"=30/SUM(I{r0+1}:I{r0+len(L)})").number_format = "0"
ws.cell(rv, 6).font = b
ws.cell(rv + 1, 1, "Clase NSR-10")
ws.cell(rv + 1, 6, f'=IF(F{rv}>760,"B",IF(F{rv}>360,"C",IF(F{rv}>=180,"D","E")))')
ws.cell(rv + 3, 1, "Fórmulas: R = √(z²+x²); tv = t·z/R; V capa = Δz/Δtv = 1/pendiente; "
        "ν = ((Vp/Vs)²−2)/(2((Vp/Vs)²−1)); Gmax = (γ/g)·Vs²; Vs30 = 30/Σ(h/Vs)")
for col, w in zip("ABCDEFGHIJKLM", [8, 9, 9, 9, 13, 13, 13, 13, 13, 14, 14, 14, 12]):
    ws.column_dimensions[col].width = w
ws.row_dimensions[5].height = 45; ws.row_dimensions[r0].height = 32
ch = ScatterChart(); ch.title = "Tiempo vertical vs profundidad"; ch.style = 13
ch.x_axis.title = "tv (ms)"; ch.y_axis.title = "z (m)"; ch.y_axis.scaling.orientation = "maxMin"
ch.x_axis.delete = False; ch.y_axis.delete = False
yz = Reference(ws, min_col=2, min_row=6, max_row=last)
for col, name in ((7, "Onda S"), (6, "Onda P")):
    s = Series(yz, Reference(ws, min_col=col, min_row=6, max_row=last), title=name)
    s.marker.symbol = "circle"; s.graphicalProperties.line.noFill = True; ch.series.append(s)
ch.height, ch.width = 12, 14
ws.add_chart(ch, "O5")
wb.save(F)
print("ok", r0, rv)
