# Sismo-geotecnia

Proyecto de recopilación, homogenización y análisis de catálogos sísmicos para Neiva, Huila. Se integran registros del Servicio Geológico Colombiano (SGC) y catálogos globales, con radios de búsqueda de 50 km y 260 km desde el Parque Central Santander.

## Estructura

```text
codigo/       Notebook de análisis y scripts auxiliares
datos/        Catálogos originales e intermedios
figuras/      Mapas y gráficos utilizados en el informe
informe/      Informe LaTeX, bibliografía y PDF
recursos/     Material de referencia
requirements.txt
```

## Entorno

Se requiere Python 3.10 o posterior y las dependencias de `requirements.txt`. Para crear un entorno virtual desde la raíz:

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install -r requirements.txt
```

Para ejecutar el análisis, abre `codigo/Codigo_madre.ipynb` en VS Code, selecciona el entorno Python y ejecuta las celdas en orden. El notebook procesa los CSV originales de SGC y USGS, combina los catálogos, calcula `escala_mw`, aplica el criterio de réplicas y presenta estadísticas y gráficos.

## Salidas del notebook

Los archivos de `datos/datos depurados/` corresponden a estas etapas:

| Archivo | Contenido |
| --- | --- |
| `Sismos_<radio>_Depurados.csv` | Catálogo combinado usado como entrada a la conversión. |
| `Sismos_<radio>_Depurados_mw.csv` | Registros con una magnitud `escala_mw` directa o calculada y metadatos de conversión. |
| `Sismos_<radio>_Depurados_onlymw.csv` | Catálogo tras el filtro de magnitud y la eliminación de réplicas. |

`<radio>` corresponde a `50KM` o `260KM`. Las estadísticas descriptivas y los gráficos del notebook se generan al ejecutar las celdas; las salidas guardadas incluyen:

- `Sismos_50KM_Depurados.csv`, `Sismos_50KM_Depurados_mw.csv` y `Sismos_50KM_Depurados_onlymw.csv`.
- `Sismos_260KM_Depurados.csv`, `Sismos_260KM_Depurados_mw.csv` y `Sismos_260KM_Depurados_onlymw.csv`.

El criterio propio de declustering usa una separación máxima de 10 km y 24 horas, dentro del mismo día. En la última ejecución guardada se reportaron 110 eventos finales en 50 km y 2.627 en 260 km; se identificaron 0 y 247 réplicas, respectivamente.

### Umbral de magnitud

La última ejecución guardada del notebook aplica `Mw >= 3`. Por eso, el archivo regional `Sismos_260KM_Depurados_onlymw.csv` incluye 33 eventos con `Mw = 3`. El informe y los gráficos `catalogo_mw_mayor_3_50km.png` y `catalogo_mw_mayor_3_260km.png` usan el subconjunto estricto `Mw > 3`: 110 eventos para 50 km y 2.594 para 260 km. Tener en cuenta esta diferencia al comparar las salidas.

Los resultados reflejan solo los registros que obtuvieron una magnitud de momento utilizable con las conversiones aplicadas. No equivalen al catálogo completo descargado ni estiman por sí solos la completitud sísmica o los parámetros Gutenberg--Richter.

## Scripts auxiliares

Los scripts de `codigo/` incluyen utilidades para descargar catálogos, consolidar fuentes y generar mapas. El flujo de análisis descrito en este README corresponde al notebook `Codigo_madre.ipynb`.

## Compilar el informe

Se requiere MiKTeX u otra distribución LaTeX, con `pdflatex` y `biber` en el `PATH`. Desde PowerShell:

```powershell
Push-Location informe
pdflatex -interaction=nonstopmode -halt-on-error main.tex
biber main
pdflatex -interaction=nonstopmode -halt-on-error main.tex
pdflatex -interaction=nonstopmode -halt-on-error main.tex
Pop-Location
```

El PDF se genera en `informe/main.pdf`. La bibliografía está en `informe/biblio.bib`.

## Archivos generados

Los archivos auxiliares de LaTeX y los registros de compilación se excluyen mediante `.gitignore`. Los catálogos, gráficos y el PDF son productos derivados de los datos y del análisis.
