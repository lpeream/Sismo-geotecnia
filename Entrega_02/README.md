# Entrega 02

Esta entrega conserva la estructura de carpetas de `Entrega_01`: `codigo`,
`datos`, `figuras`, `informe` y `recursos`, además de `requirements.txt`.

El informe está en `informe/entrega_02.tex` y su PDF compilado en
`informe/entrega_02.pdf`. Incluye un capítulo inicial de descripción del
proyecto y procedimiento, con el contexto de asociación sísmica y proximidad a
fallas de la Entrega 1. El capítulo de recurrencia ejecuta primero la
completez por intervalos y el ajuste Gutenberg--Richter en escalas local
(50 km) y regional (260 km). Incluye regresiones segmentadas, R², comparación
de tasas crudas y corregidas, gráficos por sismofuente y una curva conjunta
con los ajustes corregidos estimables y las aproximaciones de baja muestra.
No se descartan intervalos por tener menos de 10 eventos: si no se estima la
completez, se conserva su tasa cruda de referencia, identificada como no
corregida y sin incluirla en el ajuste final corregido. El notebook reproducible está en
`codigo/codigo_entrega2.ipynb`; usa la copia del catálogo en
`recursos/Qgis/datos/` y genera figuras en `figuras/` y tablas CSV en
`datos/datos estadistica/`. Entre ellas, `completez_por_tramo_magnitud.csv`,
`tabla_recurrencia_gr_por_tramo.csv` y
`conteos_incrementales_ventana_umbral.csv` conservan las métricas de los
ajustes y el desglose de conteos para cada ventana completa;
`ajustes_crudos_y_corregidos_gutenberg_richter.csv` compara los parámetros
antes y después de corregir por completez. La atenuación solo se ejecuta
después de completar la recurrencia y sus gráficas. Esta entrega
también contiene copias de todos los CSV de `Entrega_01/datos/`, con la misma
estructura de subcarpetas. El inventario general está en `datos/README.md`.

El análisis exploratorio de atenuación de Neiva está en
`codigo/analisis_atenuacion.py`. Sus entradas editables (sitio, coordenadas,
clase de sitio, PGA o `Sa(0.5 s)`, fuentes, intervalos y profundidad) están al
inicio del módulo. Genera predicciones por evento, resúmenes por fuente y
magnitud, figuras con distancia en escalas lineal y logarítmica por sismofuente
y un JSON con la configuración aplicada. La clase III
representa la hipótesis de suelo medio; esos resultados no son movimientos
en roca ni una GMPE validada localmente.

Para compilar desde la raíz del proyecto:

```powershell
Push-Location Entrega_02\informe
pdflatex -interaction=nonstopmode -halt-on-error entrega_02.tex
biber entrega_02
pdflatex -interaction=nonstopmode -halt-on-error entrega_02.tex
pdflatex -interaction=nonstopmode -halt-on-error entrega_02.tex
Pop-Location
```
