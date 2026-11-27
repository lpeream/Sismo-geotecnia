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

- `recursos/`: enunciado, datos y los dos Excel de trabajo, que son la base del
  informe (sus gráficas se exportan a `figuras/`):
  - `DATOS_DOWN_HOLE_Interpretación.xlsx`: registros del down hole y hojas
    `GRÁFICA VP` y `GRÁFICA VS` (registros apilados). En `GRÁFICA VS` están los
    tiempos leídos por canal (T23:U53), el Vs adoptado por estrato (V18:V19),
    las rectas de esas pendientes dibujadas sobre la gráfica (T56:V62) y una
    verificación por regresión (V15:V16).
  - `ensayos_analisis_sitio_Neiva (1).xlsx`: clasificación, granulometría,
    límites y ensayos dinámicos. En `Sheet1` están Vs, γ, ρ y Gmax por estrato
    y las curvas G/Gmax y ξ; en `Perfil Vs`, el perfil y el Vs30; en
    `Seed-Idriss E1` (borrador) la arena frente a Seed e Idriss (1970); en
    `Vucetic-Dobry E2` la arcilla frente a Vucetic y Dobry (1991), con los
    valores de las curvas para IP 0 a 200.
  - `DATOS DOWN HOLE 2026.xlsx` (original, no se modifica), enunciado y
    trabajo de referencia.
- `datos/ensayos_laboratorio/`: ensayos del enunciado transcritos a CSV.
- `datos/sismos/` y `datos/deepsoil/`: para los puntos 4 a 6.
- `figuras/`: figuras del informe (las `*_excel.png` salen de los Excel).
- `informe/`: informe en LaTeX (`entrega_03.tex`, se edita directamente) y su PDF.

## Notas para continuar (puntos 3 a 7)

- **Perfil adoptado (puntos 1 y 2):** dos estratos cortados en 15 m, una
  pendiente por estrato en la gráfica de la onda S. E1, arena con grava y limo
  (SP-SM) de 0 a 15 m: Vs = 150 m/s, γ = 18.6 kN/m3, ρ = 1896 kg/m3,
  Gmax = 42.7 MPa. E2, arcilla CL (CH en la base) de 15 a 30 m: Vs = 385 m/s,
  γ = 18.7 kN/m3, ρ = 1906 kg/m3, Gmax = 282.5 MPa. Vs30 = 216 m/s, perfil
  tipo D. La columna resonante queda por debajo del Gmax del down hole en los
  dos estratos (0.94 y 0.94).
- **Curvas de degradación (punto 2):** la arena (E1) queda con los ensayos
  (columna resonante de 10 m y triaxial cíclico de 5 m, de 0.002 a 10 %). La
  arcilla (E2) se completa con la curva de Vucetic y Dobry (1991) de IP = 50
  para G/Gmax y amortiguamiento (Tabla 12 del informe): es la que mejor pasa
  por los puntos de la columna resonante, aunque el IP real es 18 (efecto del
  confinamiento de 400 kPa). Valores leídos de Kramer (1996), Figuras 6.42 y
  6.45; las de IP 0, 15 y 30 de G/Gmax se verificaron con la tesis
  `recursos/TTS1127.pdf`. Las curvas llegan hasta cerca de 1 %.
- **Amortiguamiento mínimo:** la columna resonante mide valores muy bajos a
  deformaciones pequeñas (0.02 % en la arcilla, 0.25 % en la arena). La curva
  de IP = 50 arranca en cerca de 1.3 %; revisar qué valor mínimo se usa en el
  programa de respuesta de sitio (la referencia usó 5 % en todas las capas).
- **Programa:** DEEPSOIL ya no se descarga de la Universidad de Illinois; desde
  julio de 2025 es RSSeismic (Rocscience). Confirmar con el profesor.
- **Muestras M3 y M12:** clasificadas como SP-SM, como dice el profesor:
  pasa por el No. 4 el 52.9 % del total. Ojo: con la fracción gruesa (ASTM
  D2487) darían GP-GM (la grava es el 52 % de la fracción gruesa); están en
  el límite. Tienen los mismos datos.
- **Errores por corregir en el Excel de ensayos:** (1) `Sheet1` D14: el
  amortiguamiento de la columna resonante a 0.05 % es 2.05, no 1.05;
  (2) `Granulometría` N4 y N6: D10 de M3 y M12 es 0.078 mm, no 0.8;
  (3) `Resumen clasificación` K: la columna "PI" es la línea A, 0.73(LL−20);
  (4) `Triaxial cíclico M5` I: τ = G·γ/100 (31.8 kPa en el primer ensayo).
