# Entrega 02

Esta entrega conserva la estructura de carpetas de `Entrega_01`: `codigo`,
`datos`, `figuras`, `informe` y `recursos`, además de `requirements.txt`.

El informe está en `informe/entrega_02.tex` y su PDF compilado es
`informe/entrega_02.pdf`. Se organiza en cuatro secciones: sismofuentes,
recurrencia y completez, modelos de atenuación, y conclusiones para integrar
la amenaza. El Paso 2 se enfoca exclusivamente en el catálogo regional
(260 km); el radio local (50 km) se conserva como contexto espacial en el
mapa.
La completez y las tasas anuales por sismofuente regional se agrupan en
intervalos semiabiertos de magnitud 1.0, desde Mw 3.0 hasta 8.0. El script
`codigo/analisis_completez_agregada.py` reproduce la tabla regional depurada,
el gráfico pareado de conteos y tasas corregidas y la curva regional
Gutenberg--Richter por sismofuente. Para la recurrencia Gutenberg--Richter se
mantiene el ajuste por umbral en intervalos semiabiertos de 0.5, desde Mw 3.0
hasta 7.5. Para magnitudes inferiores a Mw 5.0 se estima la ventana completa
con un ajuste segmentado; para Mw >= 5.0 se adopta el periodo histórico
1917--2025. El notebook reproducible de recurrencia está en
`codigo/codigo_entrega2.ipynb`; usa la copia del catálogo en
`recursos/Qgis/datos/` y genera figuras en `figuras/` y tablas CSV en
`datos/datos estadistica/`. Entre ellas, `completez_por_tramo_magnitud.csv`,
`tabla_recurrencia_gr_por_tramo.csv` y
`conteos_incrementales_ventana_umbral.csv` conservan las métricas de los
ajustes y el desglose de conteos para cada ventana completa;
`ajustes_crudos_y_corregidos_gutenberg_richter.csv` compara los parámetros
antes y después de corregir por completez. Esta entrega también contiene
copias de los CSV de `Entrega_01/datos/`, con la misma estructura de
subcarpetas. El inventario general está en `datos/README.md`.

El análisis de atenuación de Neiva está en
`codigo/analisis_atenuacion.py`. Sus entradas editables (sitio, coordenadas,
clase de sitio, PGA o `Sa(0.5 s)`, fuentes, intervalos y profundidad) están al
inicio del módulo. Genera predicciones por evento, resúmenes por fuente y
magnitud, una figura comparativa logarítmica para las tres sismofuentes y un
JSON con la configuración aplicada. La clase III representa el escenario de
suelo medio adoptado por el modelo Zhao et al. (2006).

Para compilar desde la raíz del proyecto:

```powershell
Push-Location Entrega_02\informe
pdflatex -interaction=nonstopmode -halt-on-error entrega_02.tex
biber entrega_02
pdflatex -interaction=nonstopmode -halt-on-error entrega_02.tex
pdflatex -interaction=nonstopmode -halt-on-error entrega_02.tex
Pop-Location
```
