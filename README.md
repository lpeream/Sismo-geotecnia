# Sismo-geotecnia

Proyecto para la recopilacion y organizacion de catalogos sismicos alrededor de Neiva, Huila. Se utilizan datos del Servicio Geologico Colombiano (SGC) y del United States Geological Survey (USGS), con radios de busqueda de 50 km y 260 km.

## Estructura

```text
codigo/       Scripts para descargar, consolidar y graficar datos
datos/        Catalogos originales y catalogos consolidados
figuras/      Imagenes y mapas usados por el informe
informe/      Documento LaTeX, bibliografia y PDF generado
recursos/     Material de referencia
requirements.txt
```

## Requisitos

- Python 3.10 o posterior.
- MiKTeX u otra distribucion LaTeX.
- `pdflatex` y `biber` disponibles en el `PATH`.

### Entorno Python

Desde la raiz del proyecto, opcionalmente crea un entorno virtual:

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install -r requirements.txt
```

## Compilar el informe

El documento usa `biblatex` con backend `biber`. Desde PowerShell ejecuta:

```powershell
cd "C:\Proyectos_de_codigo\Universidad\Sismo\Sismo-geotecnia\informe"
pdflatex -interaction=nonstopmode -halt-on-error main.tex
biber main
pdflatex -interaction=nonstopmode -halt-on-error main.tex
pdflatex -interaction=nonstopmode -halt-on-error main.tex
```

El resultado queda en:

```text
informe\main.pdf
```

Para abrirlo desde PowerShell:

```powershell
Start-Process .\main.pdf
```

Si solo se modifico el texto y no la bibliografia, normalmente basta con ejecutar `pdflatex` dos veces. Cuando se agreguen o cambien referencias en `biblio.bib`, debe ejecutarse tambien `biber main`.

## Scripts

Los scripts esperan ejecutarse desde la carpeta `codigo/`.

Descargar catalogos SGC y USGS para 50 km y 260 km:

```powershell
cd codigo
python descargar_catalogos.py 50 260
```

Consolidar las fuentes en archivos CSV:

```powershell
python consolidar_fuentes.py
```

Depurar replicas y precursores con ventanas espacio-tiempo de
Gardner--Knopoff:

```powershell
python depurar_replicas.py
```

El script genera los catalogos `*-DEPURADO-GK.csv` y conserva la trazabilidad
de los eventos excluidos en archivos `*-REPLICAS-GK.csv`. Las graficas de
frecuencia--magnitud se construyen despues con:

```powershell
python analisis_frecuencia_magnitud.py
```

Generar el mapa de la zona de estudio:

```powershell
python mapa_zona_estudio.py
```

Los scripts escriben resultados en `datos/datos originales/`, `datos/datos depurados/` y `figuras/`.

## Bibliografia

Las referencias se encuentran en `informe/biblio.bib` y se citan desde `informe/main.tex` mediante `biblatex`.

## Nota sobre archivos generados

Los archivos auxiliares de LaTeX y los registros de compilacion estan excluidos mediante `.gitignore`. El codigo fuente, los datos, las figuras, la bibliografia y el PDF del informe se conservan en el repositorio.
