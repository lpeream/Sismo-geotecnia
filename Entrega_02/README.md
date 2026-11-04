# Entrega 02

Esta entrega conserva la estructura de carpetas de `Entrega_01`: `codigo`,
`datos`, `figuras`, `informe` y `recursos`, además de `requirements.txt`.

El informe está en `informe/entrega_02.tex` y su PDF compilado en
`informe/entrega_02.pdf`. Incluye un capítulo inicial de descripción del
proyecto y procedimiento, con el contexto de asociación sísmica y proximidad a
fallas de la Entrega 1. El capítulo de recurrencia presenta el análisis
Gutenberg--Richter por sismofuente y la estimación de completez por intervalos
de magnitud, incluyendo regresiones segmentadas, R², tasas normalizadas por
exposición y ajuste Gutenberg--Richter. El notebook reproducible está en
`codigo/codigo_entrega2.ipynb`; usa la copia del catálogo en
`recursos/Qgis/datos/` y genera figuras en `figuras/` y tablas CSV en
`datos/datos estadistica/`. Entre ellas, `completez_por_tramo_magnitud.csv`,
`tabla_recurrencia_gr_por_tramo.csv` y
`conteos_incrementales_ventana_umbral.csv` conservan las métricas de los
ajustes y el desglose de conteos para cada ventana completa. Esta entrega
también contiene copias de todos los CSV de `Entrega_01/datos/`, con la misma
estructura de subcarpetas. El inventario general está en `datos/README.md`.

El análisis exploratorio de atenuación de Neiva está en
`codigo/analisis_atenuacion.py`. Sus entradas editables (sitio, coordenadas,
clase de sitio, PGA o `Sa(0.5 s)`, fuentes, intervalos y profundidad) están al
inicio del módulo. Genera predicciones por evento, resúmenes por fuente y
magnitud, una figura ampliada por sismofuente y un JSON con la configuración aplicada. La clase III
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
