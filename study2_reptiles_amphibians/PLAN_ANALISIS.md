# Plan de análisis — Estabilidad proteica predicha y longevidad en reptiles y anfibios

Escrito el 9 de octubre de 2026, antes de descargar proteomas, ortólogos o estructuras y antes de cualquier prueba de
asociación. Se basa en el pipeline del estudio en mamíferos (carpeta `PROTEINAS_LONGEVIDAD`), con las lecciones de
ese trabajo incorporadas al diseño. Plan interno; para valor de prerregistro debe depositarse con su hash SHA256.

**Información vista antes de escribir el plan:** para evaluar la viabilidad se cruzó la lista de genomas anotados del
NCBI con AnAge y se vieron la vida máxima, la masa y la calidad de esos registros (48 especies con longevidad
elegible). No se ha visto ningún dato de proteínas de estas especies.

## 1. Pregunta

¿Las sustituciones de aminoácidos acumuladas desde el ancestro común de los tetrápodos hasta cada reptil o anfibio
actual tienen, en promedio, un efecto predicho distinto sobre la estabilidad de plegamiento en especies longevas
(p. ej., tortugas) que en especies de vida corta (p. ej., muchos lagartos y ranas), tras ajustar por masa y filogenia?
Prueba bilateral: hipótesis de estabilidad (menos desestabilizantes en longevas) frente a teoría casi neutral (más
desestabilizantes en especies grandes y longevas con menor tamaño poblacional efectivo).

## 2. Lecciones del estudio en mamíferos incorporadas

1. **Plantilla estructural equidistante de todas las especies analizadas:** el modelo AlphaFold del **pez cebra**
   (*Danio rerio*), que diverge de todos los tetrápodos en el mismo nodo. No se usa la estructura de ninguna especie
   analizada.
2. **Alineamientos múltiples (MAFFT)** en lugar de alineamientos por pares anclados en una especie, para reducir el
   efecto de referencia residual observado en mamíferos. Las columnas se definen por los residuos del pez cebra
   (la plantilla), que también es equidistante.
3. Referencia evolutiva: **ancestro reconstruido** (ancestro común de anfibios y reptiles), no una especie actual.
4. Controles del predictor y control positivo (tasa de sustitución) fijados de antemano.

## 3. Datos

- **Especies:** genomas de referencia del NCBI con anotación génica (NCBI Datasets; `has_annotation`, `reference_only`)
  de Testudines (8459), Crocodylia (1294634), Lepidosauria (8504) y Amphibia (8292). Aves excluidas.
- **Fenotipo:** AnAge (descarga del 8/10/2026): vida máxima > 0, calidad "acceptable" o "high", tamaño muestral distinto
  de "tiny". Coincidencia por binomio (dos primeros términos del nombre del NCBI). Sin coincidencia: sinónimos de NCBI
  Taxonomy del mismo identificador de taxón, solo binomios exactos (misma regla E1 del estudio en mamíferos).
- **Masa adulta:** AnAge; si falta, AmphiBIO v1 (`Body_mass_g`, anfibios; Oliveira et al. 2017) o ReptTraits v1-1
  (masa corporal, reptiles; Oskyrko et al. 2024); si sigue faltando, la especie se excluye. Se registra la fuente.
- **Grupos externos para la reconstrucción ancestral:** pez cebra (*Danio rerio*), pejelagarto moteado (*Lepisosteus
  oculatus*) y celacanto (*Latimeria chalumnae*). El humano se incluye en la reconstrucción como tetrápodo adicional
  (amniota), pero no en el análisis de longevidad.
- **Filogenia:** TimeTree 5, tiempos de divergencia por pares vía API (`/api/pairwise/{taxid}/{taxid}`, campo
  `adjusted_age`) para todas las especies; árbol ultramétrico por UPGMA sobre esa matriz. Se informa la desviación de
  ultrametricidad. Especies sin datos en TimeTree se excluyen (sin sustituir por congéneres).
- **Ortólogos:** archivo `gene_orthologs` del NCBI (relaciones de ortología con GeneID humano). Se exige exactamente un
  ortólogo por especie y gen.
- **Secuencias:** proteomas y GFF de cada genoma (NCBI Datasets). Por GeneID se toma la isoforma proteica más larga;
  en humano, la marcada como RefSeq Select/MANE si existe.
- **Universo de genes:** genes humanos codificantes con proteína de 50–1.000 aa y ortólogo en pez cebra cuya proteína
  es idéntica a una entrada UniProt con modelo AlphaFold DB de longitud completa. **Candidatos:** genes humanos de
  GenAge (por GeneID). **Fondo:** genes no-GenAge en orden aleatorio fijo (semilla 20261009) hasta reunir 2.000 que
  cumplen los criterios.

## 4. Alineamientos, filtros y ancestro

- MAFFT (`--auto`) por gen con todas las especies, grupos externos y humano. Columnas = posiciones del pez cebra.
- Sitios elegibles: pLDDT ≥ 70 en el modelo del pez cebra.
- Par especie–gen válido si cubre ≥ 80 % de los sitios elegibles. Gen elegible si es válido en ≥ 60 % de las especies
  analizadas y tiene al menos dos de los tres peces. Especie elegible si es válida en ≥ 70 % de los genes elegibles.
  (Filtro de genes con todas las especies, luego especies, luego genes con las especies finales.)
- IQ-TREE 2.4.0 (`-te` con el árbol de TimeTree podado, LG+G4, `-asr`): probabilidades posteriores en el ancestro común
  de anfibios y reptiles (tetrápodos).
- Máscara de errores: ventanas de 10 columnas con ≥ 5 diferencias respecto del ancestro (MAP).

## 5. Estabilidad y métricas (idénticas al estudio en mamíferos)

ThermoMPNN (pesos `thermoMPNN_default.pt`) con la secuencia ancestral enhebrada en el esqueleto del pez cebra.
ΔΔG_i(s) = f_i(b) − Σ_a P(anc_i = a)·f_i(a); K_i(s) = 1 − P(anc_i = b). Por especie: **M_s = ΣΔΔG/ΣK** (principal),
N_s = ΣΔΔG/ΣL, R_s = ΣK/ΣL.

## 6. Estadística

PGLS con λ de Pagel (MV en [0,1]) sobre el árbol de TimeTree; respuesta log10(vida máxima); predictores métrica
estandarizada + log10(masa). Bootstrap paramétrico de la razón de verosimilitudes (999; λ reestimado).

**Pruebas principales (bilaterales; Bonferroni α = 0,025):**
- **P1:** M_s con todos los genes elegibles.
- **P2:** M_s con genes GenAge; prueba competitiva frente a 1.000 conjuntos aleatorios de genes de fondo del mismo
  tamaño (p empírico bilateral) y prueba autocontenida.

**Controles del predictor:** C1 (sustituciones observadas < alternativas), C2 (enterrados > expuestos), C3 (Pro/Gly en
hélices y láminas). **Control positivo:** R_s ~ vida máxima (se espera β < 0).

**Secundarios y sensibilidad (exploratorios):** N_s; prueba casi neutral M_s ~ log10 masa; solo reptiles; solo
anfibios; sin tortugas; solo escamosos; solo masas de AnAge; sin covariable de masa; λ = 1; MCO con HC3; sin máscara;
pLDDT ≥ 90; sitios enterrados / expuestos; ajuste por log10 R_s. **Plantilla alternativa:** celacanto (también
equidistante) para 300 genes aleatorios (si hubiera < 300 genes con modelo, pejelagarto); y, solo en reptiles, la
plantilla humana (equidistante dentro de reptiles). **Temperatura:** si ReptTraits contiene temperatura corporal con
≥ 70 % de cobertura de los reptiles analizados, modelo de reptiles con esa covariable. Pruebas por gen con FDR de
Benjamini–Hochberg (exploratorias).

## 7. Interpretación fijada de antemano

- Como en mamíferos: β > 0 compatible con la teoría casi neutral; β < 0 compatible con la hipótesis de estabilidad;
  ninguna demuestra causalidad. Sin controles correctos, el resultado no es interpretable.
- La longevidad está concentrada en pocos clados (tortugas, cocodrilos); con ~40–50 especies la potencia es limitada
  y se informará el efecto detectable.
- En ectotermos la temperatura ambiental puede influir a la vez en la estabilidad proteica y en la longevidad; si no
  se puede ajustar, se declarará como confusor no controlado.

## 8. Semilla y software

Semilla 20261009. Mismo entorno de Python que el estudio en mamíferos (`PROTEINAS_LONGEVIDAD/.venv`), ThermoMPNN e
IQ-TREE 2.4.0 copiados en `tools/`, MAFFT (binario oficial para Windows).

## Enmiendas

**E1 (9/10/2026, operativa; antes de descargar proteomas, ortólogos o estructuras).** (a) Figshare ofrece ReptTraits
v1-2 (versión corregida de v1-1); se usa v1-2. Su campo de masa es "Maximum body mass (g)" (masa máxima, no típica),
al igual que `Body_mass_g` de AmphiBIO; la sensibilidad "solo masas de AnAge" ya prevista cubre esta mezcla de
definiciones. (b) ReptTraits incluye "Mean Tb" (temperatura corporal) y "Mean Annual Temperature (°C)" del área de
distribución. Además de la sensibilidad con Mean Tb prevista, se añade otra en reptiles con la temperatura media anual,
con el mismo umbral de cobertura (≥ 70 %). Para anfibios no hay variable de temperatura comparable en AmphiBIO.

**E2 (10/10/2026, POST HOC y exploratoria: escrita después de conocer todos los resultados y antes de calcular nada de
lo que sigue).** Dos resultados quedaron sin examinar. Ninguno sustituye a las pruebas principales.

*(a) Asociación nominal con la plantilla humana en reptiles.* Sobre los 300 genes y los sitios de la sensibilidad de
plantillas (`results/template_sensitivity_species_M.csv`), en los 29 reptiles y con las plantillas humana y de pez
cebra: (1) prueba de razón de verosimilitudes por bootstrap paramétrico (999 simulaciones; semillas 20261080 y
20261081); (2) modelo con covariable de clado (tortugas + cocodrilos frente a escamosos); (3) análisis dentro de
clados: solo escamosos (n = 15) y solo tortugas + cocodrilos (n = 14); (4) ajuste por log10 de la tasa de sustitución
R (todos los genes); (5) descriptivo: correlación entre z(M_humana) − z(M_pez cebra) y log10 R. Interpretación fijada:
si β con plantilla humana se acerca a 0 y pierde la significación nominal al añadir el clado o dentro de los clados,
la asociación depende del contraste entre clados y no se interpreta como efecto de la longevidad.

*(b) Exceso de genes con P < 0,05 en las pruebas por gen (207 de 1.798; 90 esperados).* Calibración por simulación
bajo la hipótesis nula: se ajusta el modelo nulo (log10 longevidad ~ log10 masa, λ por máxima verosimilitud) a las 40
especies y se simulan 1.000 respuestas (semilla 20261082). Para cada una se repiten exactamente las pruebas por gen
(mismas especies por gen, λ por máxima verosimilitud, prueba de Wald). Se informan: media e intervalo del 95 % del
número de genes con P < 0,05 bajo la nula; P calibrado del número observado, (1 + #{recuento* ≥ 207})/1.001; P
empírico por gen, (1 + #{|t*| ≥ |t|})/1.001, con el número de genes con P empírico < 0,05 y la corrección de
Benjamini–Hochberg. Variante: las mismas pruebas con λ fijado en el valor del modelo nulo. Interpretación fijada: si el
recuento observado cae dentro del intervalo nulo, el exceso se atribuye a la calibración de las pruebas por gen y a la
dependencia entre genes, no a una asociación.
