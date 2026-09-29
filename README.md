# Del ADN a la proteína

**Bioinformática · Grado en Ciencia e Ingeniería de Datos · ULPGC**  
**Autoras:** Lucía Xufang Cruz Toste y Carlota Ayala Pérez

Esta práctica recorre la replicación, la transcripción, la traducción, el *splicing* alternativo y la relación entre secuencia y estructura proteica. Las respuestas y el código ejecutable de los seis ejercicios están en [`ejercicios-01.ipynb`](ejercicios-01.ipynb). El informe de entrega está en [`Informe_Del_ADN_a_la_Proteina.pdf`](Informe_Del_ADN_a_la_Proteina.pdf) y su fuente editable en [LaTeX](Informe_Del_ADN_a_la_Proteina.tex).

## Ejercicio 1. Replicación del ADN

**Secuencia inicial:**

- Hebra 1: 5′-ATG CCG TTA GCT-3′.
- Hebra 2: 3′-TAC GGC AAT CGA-5′.

En una ronda de replicación semiconservativa se obtienen dos moléculas, cada una con una hebra antigua y una nueva:

| Molécula hija | Hebra 5′ → 3′ | Hebra 3′ → 5′ |
| --- | --- | --- |
| 1 | ATG CCG TTA GCT (original) | TAC GGC AAT CGA (nueva) |
| 2 | ATG CCG TTA GCT (nueva) | TAC GGC AAT CGA (original) |

La **helicasa** separa las hebras; la **primasa** fabrica cebadores de ARN; la **ADN polimerasa** sintetiza ADN complementario en dirección 5′ → 3′ y corrige algunos errores; y la **ligasa** une fragmentos de ADN, como los de Okazaki. Si un error de copia escapa a la reparación y queda fijado, puede transmitirse a las células hijas. Su efecto depende de la posición y del cambio concreto; una mutación también puede ser silenciosa.

La extensión con Biopython obtiene `TACGGCAATCGA` mediante `complement()` al alinear ambas hebras en sentidos opuestos, y `AGCTAACGGCAT` mediante `reverse_complement()` al escribir la hebra complementaria en dirección 5′ → 3′. El resultado coincide con la resolución manual.

## Ejercicio 2. Transcripción del ADN a ARN

Para 5′-ATG CCT GAA TGC-3′ / 3′-TAC GGA CTT ACG-5′, la **hebra molde** es la segunda. La ARN polimerasa la lee de 3′ a 5′ y sintetiza **5′-AUG CCU GAA UGC-3′**. El transcrito coincide con la hebra codificante al sustituir T por U.

El **promotor** participa en el inicio y la regulación de la transcripción; no es posible señalar sus bases exactas con estas doce bases aisladas. La **región codificante** (CDS) es la parte que se traduce a aminoácidos. Un transcrito puede contener además regiones no traducidas (UTR), de modo que no todo el ARN transcrito es CDS.

La extensión lee [`ejercicio2_real.fasta`](ejercicio2_real.fasta) con `SeqIO.read()` y produce `AUGCCUGAAUGC`. Al invertir el papel de las hebras, la complementaria inversa tomada como nueva hebra codificante produciría `GCAUUCAGGCAU`. Es un experimento de orientación: no demuestra que el gen original se transcriba en ambos sentidos.

## Ejercicio 3. Traducción del ARNm a proteína

| Codón de 5′-AUG UAU GCU UAA-3′ | Interpretación |
| --- | --- |
| AUG | Inicio y metionina (Met) |
| UAU | Tirosina (Tyr) |
| GCU | Alanina (Ala) |
| UAA | Parada; no añade aminoácido |

La proteína es **Met–Tyr–Ala**. `Bio.Seq.translate()` devuelve `MYA*`; el asterisco representa la señal de parada, no un cuarto aminoácido.

Si el AUG inicial cambiase a GUG, ese codón codificaría valina en el código genético estándar, pero en este contexto eucariota podría perderse el inicio de traducción en esa posición; podría existir otro inicio aguas abajo. Si UAA desapareciese, la traducción podría continuar en el mismo marco hasta un codón de parada posterior, con posibles efectos en longitud y función. No se puede predecir el resultado final sin la secuencia posterior.

## Ejercicio 4. *Splicing* alternativo

En un gen hipotético con exones 1–2–3–4–5, dos combinaciones posibles son **1–2–4–5** y **1–3–5**, si se conservan señales de empalme compatibles. La inclusión de exones distintos puede alterar la secuencia y los dominios de las proteínas resultantes; también puede cambiar solo regiones no traducidas o generar ARN que no llegue a producir una proteína estable. Por eso hace falta conocer la secuencia y el marco de lectura para predecir una función concreta.

**Ejemplo real: FGFR2.** La selección de las regiones IIIb y IIIc modifica parte del dominio extracelular y puede cambiar la unión a ligandos. El notebook consulta los transcritos de FGFR2 anotados en Ensembl y **compara dos variantes concretas** de RefSeq, `NM_022970` (IIIb) y `NM_000141` (IIIc). Descarga sus registros GenBank, extrae las CDS y muestra las longitudes de ARNm y proteína y los tramos proteicos distintos. Las versiones y longitudes se imprimen durante la ejecución: dependen de la anotación recuperada. El número total de transcritos de Ensembl no equivale al número demostrado de proteínas funcionales, y la comparación de secuencias por sí sola no atribuye todas las diferencias a un único exón.

El *splicing* aumenta la variedad de productos de ARN y, en algunos casos, de proteínas sin añadir genes. Su regulación puede variar entre tejidos.

## Ejercicio 5. Secuencia y estructura de las proteínas

En **Met–Ile–Ser–Gly–Val–Lys–His**, el extremo **N-terminal** es Met y el **C-terminal** es His. El orden de los residuos condiciona las interacciones que permiten el plegamiento y, por tanto, la función; el entorno celular también importa.

Si un residuo hidrofóbico enterrado en el núcleo se sustituyera por otro hidrofílico, podría disminuir la estabilidad o cambiar el plegamiento. El efecto exacto dependería del residuo y de su entorno: no toda sustitución causa pérdida de función.

Como extensión se representa la **GFP (PDB `1EMA`)** con `py3Dmol`. En ella se observa un barril de láminas beta que rodea la región central del cromóforo. Una mutación que alterase interacciones necesarias para mantener esa estructura podría afectar al plegamiento o a la fluorescencia. La visualización interactiva se abre desde el notebook y usa primero el archivo local [`pdb1ema.ent`](pdb1ema.ent).

## Ejercicio 6. Actividad integradora

Se emplea el registro humano de insulina **NCBI RefSeq `NM_000207.3`**. El notebook descarga su registro GenBank, extrae la **CDS de 333 nucleótidos** y la guarda en `insulina_CDS_NM_000207_3.fasta` al ejecutarse. Esta CDS del transcrito ya procesado sirve como modelo compacto para los tres pasos: no representa un locus genómico con intrones ni una simulación molecular completa de una horquilla de replicación.

1. **Replicación:** se calcula la hebra complementaria y se muestran dos moléculas hijas con una hebra original y otra nueva.
2. **Transcripción:** la hebra molde 3′ → 5′ determina el ARNm 5′ → 3′, equivalente a transcribir la hebra codificante mediante Biopython.
3. **Traducción:** el ARNm produce una **preproinsulina de 110 aminoácidos**, antes del procesamiento que da lugar a la insulina madura. El codón de parada no se incluye como aminoácido.

El pipeline muestra la longitud completa procesada y abrevia solo la **presentación** de las secuencias largas. En el registro de ejecución anterior, la salida comenzaba con ADN `ATGGCCCTGTGGATG...`, ARNm `AUGGCCCUGUGGAUG...` y proteína `MALWMRLLPLLALLALW...`.

**Reflexión:** un error de replicación que se fije en el ADN puede persistir en una línea celular y afectar a futuros transcritos; un error aislado de transcripción o traducción suele afectar a moléculas concretas, sin modificar el ADN. Esto no implica que toda mutación cambie una proteína ni que toda proteína derivada de un alelo mutado sea necesariamente defectuosa.

## Archivos y ejecución

| Archivo | Contenido |
| --- | --- |
| [`ejercicios-01.ipynb`](ejercicios-01.ipynb) | Respuestas y código de los seis ejercicios; ejecutar para renovar las consultas |
| [`ejercicio2_real.fasta`](ejercicio2_real.fasta) | Secuencia de ejemplo para la transcripción |
| [`pdb1ema.ent`](pdb1ema.ent) | Estructura PDB de la GFP utilizada en la visualización |
| [`Informe_Del_ADN_a_la_Proteina.pdf`](Informe_Del_ADN_a_la_Proteina.pdf) | Informe académico de cuatro páginas como máximo |
| [`Informe_Del_ADN_a_la_Proteina.tex`](Informe_Del_ADN_a_la_Proteina.tex) | Fuente LaTeX editable del informe |
| `insulina_CDS_NM_000207_3.fasta` | Se genera al ejecutar el ejercicio 6 desde RefSeq |
| [`requirements.txt`](requirements.txt) | Dependencias de Python |
| [`Presentacion.md`](Presentacion.md) | Apoyo visual y guion para la exposición |
| `README.md` | Respuestas, métodos y guía de uso |

Se necesita Python 3 y las bibliotecas `biopython`, `requests`, `py3Dmol`, `jupyter`. En la carpeta del proyecto:

```bash
python -m pip install -r requirements.txt
jupyter notebook ejercicios-01.ipynb
```

Ejecuta las celdas en orden. Las consultas de FGFR2 a Ensembl y de insulina a NCBI requieren conexión; las anotaciones pueden cambiar. El notebook escribe `ejercicio2_real.fasta` en el directorio de ejecución, y la descarga PDB puede volver a generar `pdb1ema.ent`. La visualización 3D se aprecia mejor en un entorno Jupyter con soporte HTML.

## Fuentes

- Alberts, B. y colaboradores. *Molecular Biology of the Cell*. Garland Science, 2014.
- Ensembl, [FGFR2 humano](https://www.ensembl.org/Homo_sapiens/Gene/Summary?g=ENSG00000066468) y [servicio REST](https://rest.ensembl.org/).
- NCBI Nucleotide, [RefSeq `NM_000207.3`](https://www.ncbi.nlm.nih.gov/nuccore/NM_000207.3).
- NCBI Gene, [FGFR2](https://www.ncbi.nlm.nih.gov/gene/2263); variantes [IIIb, `NM_022970`](https://www.ncbi.nlm.nih.gov/nuccore/NM_022970) y [IIIc, `NM_000141`](https://www.ncbi.nlm.nih.gov/nuccore/NM_000141).
- RCSB PDB, [estructura `1EMA`](https://www.rcsb.org/structure/1EMA).

**Presentación oral:** preparar apoyo visual con una estructura, un esquema del proceso y la reflexión crítica. Añadir la URL real del repositorio al compartir la entrega; no se proporciona una URL en los materiales recibidos.
