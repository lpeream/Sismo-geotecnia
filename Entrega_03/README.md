# Entrega 03: análisis de sitio (20 %)

Perfil de suelo de 30 m en Neiva, apoyado sobre roca con el movimiento
determinado en la Entrega 02 (PGA de 475 años en roca = 0.196 g; espectro de
aceleración en roca con meseta de 0.481 g). El enunciado y los datos
crudos del down hole están en `recursos/`.

## Puntos del enunciado

1. Perfil estratigráfico, clasificación y propiedades dinámicas de cada estrato.
2. Curvas normalizadas G/Gmax y amortiguamiento por estrato.
3. Desplazamientos y G contra profundidad con el análisis modal de Zeevaert.
4. Tres sismos (fuente cercana, intermedia y lejana), con justificación.
5. Análisis en DEEPSOIL y comparación con Zeevaert.
6. Espectro de aceleración en superficie (DEEPSOIL) frente al de la NSR-10.
7. Potencial de licuación (Zeevaert y DEEPSOIL), con nivel freático en superficie.

## Estructura

- `codigo/`: scripts de cálculo. Orden de ejecución:
  1. `procesamiento_downhole.py`: llegadas P y S del down hole crudo, perfiles de
     Vp y Vs (canales cada 1 m, fuente a 3 m) y comparación con la referencia.
  2. `perfil_estratigrafico.py`: estratos, Vs30, Gmax y figura del perfil.
  3. `curvas_degradacion.py`: curvas G/Gmax y D por estrato.
- `datos/ensayos_laboratorio/`: ensayos del enunciado transcritos a CSV.
  - `clasificacion_muestras.csv` y `granulometrias.csv`: humedad, granulometría,
    límites y USCS informado (sondeo P-1).
  - `propiedades_muestras_dinamicas.csv`, `triaxial_ciclico_5m.csv` y
    `columna_resonante.csv`: ensayos dinámicos a 5, 10 y 20 m.
  - `triaxial_cu.csv` y `triaxial_cu_parametros.csv`: triaxial CU a 4 y 28 m.
  - `vs_downhole_leido_de_grafica.csv`: Vs del down hole **leído de la
    gráfica** del enunciado (valores aproximados, por verificar).
- `datos/downhole/`: copia de trabajo del Excel del down hole (la que leen los
  scripts; el original de `recursos/` no se modifica).
- `datos/sismos/`: registros de los tres sismos elegidos.
- `datos/deepsoil/`: archivos de entrada y resultados de DEEPSOIL.
- `datos/resultados/`: tablas generadas por los scripts.
- `figuras/`: figuras del informe.
- `informe/`: informe LaTeX y PDF.
- `recursos/`: enunciado (`TRABAJO ANALISIS LOCAL 2026-2.pdf`) y registros
  crudos del down hole (`DATOS DOWN HOLE 2026.xlsx`: hojas `VP`, 40 canales Z,
  y `VS`, 40 canales X/Y; 4000 muestras con dt = 0.5 ms; fuente a 3 m del
  pozo según el esquema de la hoja `REFERENCIAS PARA ENSAYO`).

## Notas para continuar (puntos 3 a 7)

- **Perfil adoptado (puntos 1 y 2):** dos estratos. E1, arena SP-SM de 0 a
  16.5 m: Vs = 151 m/s (tres rectas: 145, 148 y 160 m/s), gamma = 18.6 kN/m3,
  Gmax = 43.1 MPa. E2, arcilla CL (CH en la base) de 16.5 a 30 m: Vs = 380 m/s
  (una recta, canales 17 a 30), gamma = 18.7 kN/m3, Gmax = 275.7 MPa.
  Vs30 = 207 m/s, perfil tipo D. Valores en `datos/resultados/`.
- **Curva de la arcilla (para DEEPSOIL):** solo hay datos hasta gamma = 0.05 %
  (columna resonante de 20 m). El resto se completó por ahora con la forma de
  Darendeli ajustada a esos puntos (gamma_r = 0.65 %, a = 0.919). La idea es
  reemplazarla, como en el trabajo de referencia, por la curva de Vucetic y
  Dobry (1991) que trae DEEPSOIL: probar IP = 15 e IP = 30 (la arcilla tiene
  IP de 18), compararlas con los 5 puntos de la columna resonante y escoger la
  que mejor se ajuste. Luego actualizar la figura de curvas, la tabla de
  parámetros y el texto del punto 2 (está la nota en rojo).
- **Dmin de la arcilla:** el ajuste da 0 % (el ensayo mide 0.02 %). Hay que
  fijar un valor razonable para DEEPSOIL (usual 0.5 a 2 %; la referencia usó
  5 % en todas las capas).
- **Down hole:** el canal 18 de la onda S se ajustó a mano de 129 a 127 ms
  (la correlación saltaba un ciclo). Detalles en
  `figuras/lectura_downhole_S/canales_1_a_30.pdf`.
- **Muestras M3 y M12:** tienen los mismos datos; por granulometría darían
  GP-GM. Se dejaron como SP-SM; consultar con el profesor.
