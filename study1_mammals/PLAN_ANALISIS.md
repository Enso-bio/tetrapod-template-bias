# Plan de análisis — Estabilidad proteica predicha y longevidad en mamíferos placentarios

Escrito el 8 de octubre de 2026, antes de descargar fenotipos, ortólogos o estructuras y antes de cualquier
prueba de asociación. Plan interno, no un prerregistro público: para que tenga valor de prerregistro hay que
depositarlo (p. ej. en OSF) con su hash SHA256 antes de ejecutar `08_stats.py`. Cualquier cambio posterior se
anota en la sección "Enmiendas" con fecha y motivo, sin borrar el texto original.

## 1. Pregunta

¿Las sustituciones de aminoácidos acumuladas desde el ancestro placentario hasta cada especie actual tienen,
en promedio, un efecto predicho distinto sobre la estabilidad de plegamiento en especies longevas que en especies
de vida corta, después de ajustar por masa corporal y filogenia?

Dos predicciones compiten y por eso la prueba es bilateral:

- **Hipótesis de estabilidad**: las especies longevas acumulan sustituciones menos desestabilizantes (o más
  estabilizantes).
- **Teoría casi neutral**: las especies grandes y longevas tienen menor tamaño poblacional efectivo, la selección
  purificadora es menos eficiente y acumulan sustituciones *más* desestabilizantes.

## 2. Diferencias deliberadas con el estudio anterior (geometría de AlphaFold)

1. No se compara geometría de modelos por especie: las diferencias entre ortólogos (≈0,1–0,3 Å) estaban por debajo
   del error de los modelos. Se predice el efecto energético de cada sustitución sobre una única estructura.
2. La referencia no es el humano sino el **ancestro placentario reconstruido**. El humano es una especie más.
3. La métrica principal es por sustitución (calidad), no la suma (cantidad), para no confundir con la tasa de
   evolución, que depende del tiempo generacional.
4. Panel amplio con fondo aleatorio: genes de GenAge frente a genes de fondo muestreados al azar.
5. No se exige identidad UniProt–Ensembl para especies no humanas: se usan directamente los ortólogos 1:1 de Ensembl.

## 3. Datos

- **Fenotipo**: AnAge (descarga actual de genomics.senescence.info). Criterios: clase Mammalia; vida máxima y masa
  adulta > 0; calidad "acceptable" o "high"; tamaño muestral distinto de "tiny".
- **Filogenia**: árbol MCC con datación por nodos de Upham et al. (2019), mismo commit que el estudio anterior
  (740b70189c77fc168be17bd6ade204c755bcb19b); 20 árboles posteriores equiespaciados de la muestra de 100 para
  sensibilidad.
- **Especies**: genomas de mamíferos en Ensembl Compara (API REST, versión vigente) cuyo binomio (dos primeros
  términos del nombre de producción) coincide exactamente con Género_especie de AnAge y con una punta del árbol.
  Si varios genomas comparten binomio, se usa el que coincide exactamente; si ninguno, el primero en orden
  alfabético. Sin sinónimos supuestos ni transferencias (p. ej., perro→lobo).
- **Universo de genes**: proteínas humanas UniProt revisadas (Swiss-Prot) con secuencia canónica idéntica a la
  proteína canónica de Ensembl del gen, longitud 50–1.000 aa y modelo AlphaFold DB de cadena completa.
- **Conjunto candidato**: genes humanos de GenAge (descarga actual) presentes en el universo.
- **Fondo**: 2.000 genes del universo no incluidos en GenAge, muestreados al azar con semilla 20261008.
- Las 12 familias del estudio anterior (SIRT6, SOD1, SOD2, CAT, GPX1, PRDX1, HSPA8, OGG1, PARK7, ACTB, GAPDH,
  TPI1) se añaden si están en el universo; se analizan solo de forma exploratoria.

## 4. Alineamientos y ancestro

- Ortólogos 1:1 de Ensembl (homología pairwise con humano, secuencias de proteína alineadas). Se construye un
  alineamiento anclado en coordenadas humanas: columnas = posiciones humanas; inserciones respecto al humano se
  descartan.
- Pares especie–gen válidos si cubren ≥ 80 % de los sitios elegibles del gen.
- Gen elegible si tiene ortólogo válido en ≥ 60 % de las especies placentarias elegibles y al menos una especie no
  placentaria (marsupial o monotrema) como grupo externo.
- Especie elegible si tiene ortólogo válido en ≥ 70 % de los genes elegibles.
- **Reconstrucción ancestral**: IQ-TREE (`-te` con la topología de Upham podada, longitudes de rama reoptimizadas,
  modelo LG+G4, `-asr`). Estado ancestral en el nodo MRCA de los placentarios, con probabilidades posteriores
  marginales por sitio.
- **Filtro de errores de anotación/alineamiento** (fijado de antemano): en cada par especie–gen se enmascaran las
  posiciones situadas en ventanas de 10 posiciones humanas con ≥ 5 diferencias respecto del ancestro (estado MAP).

## 5. Predicción de estabilidad

- **Estructura**: modelo AlphaFold DB humano (solo esqueleto). Sitios elegibles: pLDDT ≥ 70.
- **Predictor**: ThermoMPNN (Dieckhaus et al. 2024, PNAS; pesos `thermoMPNN_default.pt`, licencia MIT). Se enhebra la
  secuencia ancestral (MAP) sobre el esqueleto humano; una pasada da, para cada sitio i, un vector de 20 valores
  e_i(a). El efecto de cambiar a→b en el contexto ancestral es ΔΔG_i(a→b) = e_i(b) − e_i(a); positivo =
  desestabilizante (convención de ThermoMPNN, se verificará con controles de cordura).
- Para la especie s en el sitio i: ΔΔG_i(s) = Σ_a P(anc_i = a) · [e_i(s_i) − e_i(a)].
- Selenocisteína (U) se trata como C; X y huecos se excluyen.

## 6. Métricas

Para cada especie s y gen g (sitios elegibles no enmascarados con residuo presente):

- D_sg = Σ_i ΔΔG_i(s)  (cambio neto predicho desde el ancestro)
- K_sg = Σ_i [1 − P(anc_i = s_i)]  (sustituciones esperadas)
- L_sg = número de sitios evaluados

Agregados por especie sobre un conjunto de genes G:

- **M_s(G) = Σ_g D_sg / Σ_g K_sg** — ΔΔG medio por sustitución. **Métrica principal.**
- N_s(G) = Σ_g D_sg / Σ_g L_sg — carga neta por sitio (secundaria; mezcla tasa y efecto).
- R_s(G) = Σ_g K_sg / Σ_g L_sg — tasa de sustitución (descriptiva y covariable de sensibilidad).

## 7. Modelo estadístico

- Unidad: especie placentaria. Respuesta: log10(vida máxima). Predictores: métrica estandarizada + log10(masa adulta).
- PGLS con λ de Pagel por máxima verosimilitud en [0,1] sobre el árbol MCC.
- Inferencia: razón de verosimilitudes con bootstrap paramétrico nulo (999 simulaciones, λ reestimado en cada una;
  p = (excedencias+1)/1000). IC95 de Wald condicionados a λ, aproximados.

### Pruebas principales (bilaterales; Bonferroni, α = 0,025 cada una)

- **P1 (global)**: M_s(todos los genes elegibles) ~ log10 vida máxima.
- **P2 (competitiva)**: M_s(GenAge). Estadístico: β estandarizado de M_s(GenAge). Nulo: distribución de β para 1.000
  conjuntos aleatorios de genes de fondo del mismo tamaño (semilla fija). p empírico bilateral
  = (nº |β_nulo| ≥ |β_obs| + 1)/1001. Se informa también la prueba autocontenida (PGLS de M_s(GenAge) solo).

### Controles de cordura (deben salir como se espera para confiar en el predictor)

- C1: ΔΔG medio de las sustituciones observadas < ΔΔG medio de las 19 alternativas posibles en los mismos sitios
  (selección purificadora).
- C2: |ΔΔG| de las alternativas mayor en sitios enterrados (RSA < 0,25 en el modelo humano) que en expuestos.
- C3: sustituciones a prolina o glicina en hélices/hojas predichas más desestabilizantes que el promedio.

### Secundarias y sensibilidad (exploratorias; se informan todas)

- N_s y R_s con el mismo modelo.
- Prueba teórica casi neutral: M_s ~ log10(masa) (sin vida máxima).
- λ = 1 fijo; MCO con HC3; 20 árboles posteriores; sin filtro de ventanas; solo pLDDT ≥ 90; solo sitios enterrados /
  solo expuestos; excluir primates; excluir humano; excluir murciélagos; solo registros AnAge de calidad "high";
  añadir log10 R_s como covariable; usar LQ (residuo alométrico) en lugar de vida máxima + masa.
- Plantilla alternativa: esqueleto del modelo AlphaFold de ratón para un subconjunto aleatorio de 200 genes;
  se informa la correlación de M_s entre plantillas.
- Predictor alternativo (si se logra ejecutar): ESM-2 650M con marginales enmascaradas en el contexto ancestral.
- Análisis por gen (D_sg/K_sg por gen, PGLS por gen, FDR de Benjamini-Hochberg sobre todos los genes): exploratorio.
- Las 12 familias del estudio anterior: PGLS por gen, exploratorio.
- Tres grupos por terciles de LQ (vida máxima / [4,88·masa^0,153], de Magalhães et al. 2007): solo descriptivo.

## 8. Interpretación fijada de antemano

- Asociación positiva de M_s con vida máxima (más desestabilizante en longevas) es compatible con la teoría casi
  neutral, no con la hipótesis de estabilidad.
- Asociación negativa es compatible con la hipótesis de estabilidad, pero no demuestra causalidad.
- Ausencia de asociación con controles de cordura correctos y potencia razonable restringe efectos grandes; sin
  controles correctos, el resultado no es interpretable.
- ThermoMPNN predice ΔΔG de mutaciones puntuales con error considerable (r ≈ 0,5–0,7 en bancos de prueba); las
  conclusiones son sobre estabilidad *predicha*.

## 9. Semilla y software

Semilla 20261008. Python 3.12, PyTorch con CUDA, ThermoMPNN (commit de la rama main descargado el 8/10/2026),
IQ-TREE (binario oficial), Biopython, statsmodels. Versiones exactas en `runtime_versions.json`.

## Enmiendas

**E1 (8/10/2026, tras descargar AnAge, árboles y lista de genomas de Ensembl; antes de descargar ortólogos y de
cualquier prueba de asociación).** Varias especies de Ensembl no coinciden con AnAge o con el árbol solo por
nomenclatura (p. ej., *Nannospalax galili* / *Spalax galili*, *Physeter catodon* / *P. macrocephalus*). Se añade una
reconciliación objetiva: para cada genoma de Ensembl sin coincidencia, se consultan los nombres del mismo identificador
de taxón en NCBI Taxonomy (nombre científico, sinónimos, nombres equivalentes y sinónimos de GenBank). Se acepta un
nombre alternativo solo si es un binomio exacto (dos términos; no se truncan trinomios) y coincide exactamente con
AnAge y/o con una punta del árbol. No se usan juicios manuales ni nombres de especies hermanas. Las asignaciones
quedan en `data/species_synonyms_ncbi.csv`.

**E2 (aclaración operativa, misma fecha).** El universo exige identidad entre la proteína canónica de UniProt y la de
Ensembl, que solo se comprueba al consultar Ensembl. Por eso el fondo se obtiene recorriendo una permutación aleatoria
fija (semilla 20261008) de los genes candidatos no-GenAge hasta reunir 2.000 que cumplen los criterios.

**E3 (8/10/2026, antes de obtener ortólogos y de cualquier prueba).** La API REST de Ensembl respondía en ~14–140 s
por gen, inviable para ~2.300 genes. Se usan los archivos masivos del FTP de Ensembl 116 (homologías de Compara y
FASTA de proteínas). Las llamadas de ortología 1:1 son las mismas. Como el FTP no trae los alineamientos por pares,
cada ortólogo se alinea con la proteína humana mediante alineamiento global con extremos libres (BLOSUM62, apertura
−10, extensión −0,5; Biopython PairwiseAligner), anclado en las posiciones humanas. El gen humano se identifica por
la proteína canónica que usa Compara, y su secuencia debe ser idéntica a la canónica de UniProt.

**E4 (8/10/2026, POST HOC: añadida después de ver P1, P2 y las sensibilidades; no sustituye a ninguna prueba
principal).** Las sensibilidades previstas mostraron que M_s con plantilla humana es mínimo en el humano y crece con
la distancia filogenética al humano (primates muy por debajo del resto), y que la plantilla de ratón desplaza a
primates y roedores en sentidos opuestos: un sesgo de compatibilidad entre el esqueleto de la plantilla y las
secuencias de especies cercanas a ella. Como análisis exploratorio se repite el cálculo con una plantilla
equidistante de todos los placentarios: el modelo AlphaFold DB de la zarigüeya (*Monodelphis domestica*, ortólogo
1:1 de Ensembl con secuencia idéntica a una entrada UniProt), con la secuencia ancestral enhebrada en las posiciones
alineadas con el humano y sitios con pLDDT ≥ 70 en ambos modelos. Se informan P1 y P2 con esta métrica, con los
mismos métodos, etiquetados como post hoc.

**E5 (9/10/2026, CONFIRMATORIA: escrita antes de descargar ningún dato de plantilla del ornitorrinco y antes de
calcular nada con él).** Objetivo: confirmar con una segunda plantilla externa, independiente de la zarigüeya, el
sesgo de plantilla descubierto en E4 y la ausencia de asociación. Plantilla: modelo AlphaFold DB del ortólogo 1:1 de
Ensembl del ornitorrinco (*Ornithorhynchus anatinus*, taxón NCBI 9258) con secuencia idéntica a una entrada UniProt;
el ornitorrinco es equidistante de todos los placentarios y no es pariente cercano de la zarigüeya (monotrema frente a
marsupial). Mismo procedimiento que E4 (secuencia ancestral enhebrada, pLDDT ≥ 70, mismo filtro de ventanas).

*Conjunto de comparación (fijado):* genes con modelo válido para ambas plantillas externas (zarigüeya y ornitorrinco);
sitios con pLDDT ≥ 70 en el modelo humano y en ambos modelos externos, alineados en los tres y no enmascarados. Sobre
esos mismos genes y sitios se calcula M_s con las tres plantillas (humana, zarigüeya, ornitorrinco). Si quedaran menos
de 500 genes comunes, se usarían los genes del ornitorrinco con sus propios sitios (no se espera).

*Definiciones:* desplazamiento de primates Δ_T = (media de M_s en primates − media en los demás placentarios) / DE de
M_s en las 57 especies, para cada plantilla T. Intervalos por bootstrap de genes: 1.000 remuestreos con reemplazo de
genes (semilla 20261059), recalculando M_s de las tres plantillas en cada remuestreo; IC95 por percentiles.

*Hipótesis confirmatorias y criterios de decisión:*
- **H5a (sesgo de plantilla):** Δ_humana − Δ_ornitorrinco < 0 (los primates parecen más estables con la plantilla
  humana). Se confirma si el límite superior del IC95 bootstrap es < 0.
- **H5b (las dos plantillas externas coinciden entre sí más que con la humana):** r(M_zarigüeya, M_ornitorrinco) −
  r(M_humana, M_ornitorrinco) > 0 (correlaciones de Pearson entre especies). Se confirma si el límite inferior del IC95
  bootstrap es > 0.
- **H5c (el nulo se replica):** PGLS de P1 con M_s de la plantilla del ornitorrinco (mismo modelo, 999 simulaciones,
  semilla 20261060). Se considera replicado si p bootstrap ≥ 0,025 y el IC95 de β incluye 0.

*Se informan además (descriptivos):* el rango del humano entre las 57 especies con cada plantilla, la prueba
competitiva de GenAge con la plantilla del ornitorrinco (1.000 conjuntos, semilla 20261061) y β de P1 con las tres
plantillas sobre los sitios comunes. No se añadirá ningún otro análisis a esta enmienda tras ver los resultados.

**E6 (10/10/2026, POST HOC y exploratoria: escrita después de conocer E5 y antes de calcular nada de lo que sigue).**
En E5 persistió con la plantilla del ornitorrinco un desplazamiento de primates (Δ = −0,62 DE) y el humano siguió
siendo la especie con M_s más bajo. En E4/E5 los alineamientos estaban anclados en las posiciones humanas y los sitios
exigían pLDDT ≥ 70 en el modelo humano. Objetivo: comprobar si esa ventaja residual desaparece cuando ningún paso usa
la secuencia ni el modelo humanos como referencia.

*Procedimiento (fijado):* genes = los 1.341 genes comunes de E5; especies = las 57 placentarias y los 6 grupos externos
del análisis principal. Para cada gen se reúnen las proteínas de Compara 1:1 (FTP de Ensembl 116, ya descargado) de
todas las especies, incluida la humana como una más; la del ornitorrinco es la idéntica al modelo AlphaFold de E5. Se
alinean con MAFFT 7.526 (`--auto`, un hilo) y se conservan las columnas con residuo del ornitorrinco. Sitios
elegibles: pLDDT ≥ 70 en el modelo del ornitorrinco. Par especie–gen válido si cubre ≥ 80 % de los sitios elegibles;
gen retenido si es válido en ≥ 60 % de las 57 placentarias; no se refiltran especies. Ancestro placentario con IQ-TREE
2.4.0 (`-te`, árbol MCC podado a las especies válidas, `-m LG+G4 -asr -nt 1 -seed 20261008`; nodo = MRCA de los
placentarios). ThermoMPNN con el ancestro MAP enhebrado en el modelo del ornitorrinco, misma máscara de ventanas y
M_s = ΣΔΔG/ΣK.

*Resultados que se informan:* Δ con esta métrica (IC95 por bootstrap de 1.000 remuestreos de genes, semilla
20261070); diferencias pareadas frente a Δ de E5 con plantilla humana y con plantilla de ornitorrinco anclada en el
humano, en los genes retenidos; rango del humano entre las 57 especies; correlación con M_s de E5; P1 con esta métrica
(PGLS, 999 simulaciones, semilla 20261071). Variante de descomposición: la misma métrica restringida a las columnas
cuyo residuo humano alineado tiene pLDDT ≥ 70 en el modelo humano. *Interpretación fijada:* si el IC95 de Δ incluye 0,
la ventaja residual se atribuye al anclaje en el humano (alineamiento o definición de sitios, según la variante); si
Δ < 0 con IC95 que excluye 0, la ventaja persiste sin anclaje humano y su origen sigue abierto (diferencia real de los
primates u otro sesgo, p. ej. de tasa). Todo se etiqueta como post hoc.
