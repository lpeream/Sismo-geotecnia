# Sismo-geotecnia

Proyecto de recopilación, homogenización y análisis de catálogos sísmicos para Neiva, Huila. Se integran registros del Servicio Geológico Colombiano (SGC) y del catálogo ComCat del USGS, con radios de búsqueda de 50 km y 260 km desde el Parque Central Santander. El flujo produce mapas, estadísticas descriptivas y ajustes de recurrencia Gutenberg--Richter. El informe final es `informe/entrega_01.pdf`.

## Estructura

```text
codigo/       Notebook de análisis y scripts auxiliares
datos/        Catálogos originales e intermedios
figuras/      Mapas y gráficos utilizados en el informe
informe/      Fuente LaTeX, bibliografía y PDF de entrega
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

La última ejecución guardada del notebook aplica `Mw >= 3`: 110 eventos finales para 50 km y 2.627 para 260 km. El catálogo regional incluye 33 eventos con `Mw = 3`. El ajuste Gutenberg--Richter usa `Mc = 3.05`; por eso, su muestra es de 110 eventos para 50 km y 2.594 para 260 km. Los parámetros calculados fueron `b = 0.666`, `a = 4.117`, `R2 = 0.958` para 50 km y `b = 0.728`, `a = 5.813`, `R2 = 0.982` para 260 km. Los eventos regionales de `Mw = 3` se conservan en el catálogo final, pero quedan bajo `Mc` y no entran en ese ajuste.

Los resultados reflejan solo los registros que obtuvieron una magnitud de momento utilizable con las conversiones aplicadas. No equivalen al catálogo completo descargado. Los ajustes de recurrencia son preliminares y deben evaluarse junto con la completitud, la independencia de los eventos y la sensibilidad a los parámetros empleados.

## Scripts auxiliares

Los scripts de `codigo/` incluyen utilidades para descargar catálogos, consolidar fuentes y generar mapas. El flujo de análisis descrito en este README corresponde al notebook `Codigo_madre.ipynb`.

## Compilar el informe

Para compilar el informe no se necesita ejecutar Python ni el notebook. Se requiere una distribución LaTeX, como MiKTeX o TeX Live, con `pdflatex` y `biber` disponibles en el `PATH`. El proyecto usa la clase `report` y los paquetes LaTeX `geometry`, `babel` (español), `amsmath`, `amssymb`, `graphicx`, `booktabs`, `float`, `caption`, `subcaption`, `fancyhdr`, `csquotes`, `hyperref` y `biblatex` con estilo APA. MiKTeX puede instalar automáticamente los paquetes faltantes; en TeX Live deben estar incluidos en la instalación.

El fuente y los recursos requeridos son:

- `informe/entrega_01.tex`: documento LaTeX.
- `informe/biblio.bib`: referencias procesadas por `biber`.
- `figuras/Logo_1.png` y `figuras/Logo_encabezado.png`: logotipos de portada y encabezado.
- `figuras/Mapa_Sismos_No_Depurados.png` y `figuras/Mapa_Sismos_Depurados.png`: mapas comparativos incluidos horizontalmente.
- `figuras/gutenberg_richter_50km_onlymw.png` y `figuras/gutenberg_richter_260km_onlymw.png`: gráficos de recurrencia.

Las rutas de figuras dentro del `.tex` son relativas a `informe/`; conserva esta estructura de carpetas. Desde PowerShell, en la raíz del repositorio, ejecuta:

```powershell
Push-Location informe
pdflatex -interaction=nonstopmode -halt-on-error entrega_01.tex
biber entrega_01
pdflatex -interaction=nonstopmode -halt-on-error entrega_01.tex
pdflatex -interaction=nonstopmode -halt-on-error entrega_01.tex
Pop-Location
```

La secuencia ejecuta LaTeX antes y después de `biber` para actualizar las citas, referencias cruzadas, tabla de contenido e índice de figuras. El PDF se genera en `informe/entrega_01.pdf`.

Para volver a generar los mapas o gráficos a partir de los catálogos sí se necesita Python, las dependencias de `requirements.txt` y ejecutar las celdas correspondientes de `codigo/Codigo_madre.ipynb`. La compilación LaTeX usa las imágenes existentes y no descarga mapas base.

## Archivos generados

Los archivos auxiliares de LaTeX (`.aux`, `.bbl`, `.bcf`, `.blg`, `.lof`, `.log`, `.out`, `.run.xml` y `.toc`) se generan durante la compilación y se excluyen mediante `.gitignore`. Los catálogos, gráficos y el PDF son productos derivados de los datos y del análisis.
