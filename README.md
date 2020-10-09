The structural variation pipeline calls structural variants in Illumina
paired-end reads from whole genome mouse data relative to _GRCm38 (mm10)_. The
pipeline maps reads to the reference, and then passes mapped reads to four
structural variant callers: `Breakdancer`, `Lumpy`, `Delly`, and `Manta`. Structural
variant calls from the four individual callers are then merged with `Survivor`.

Calls are merged by type, within a +/- 1000bp buffer around each call. Merged
calls provide the number of callers supporting each call.

Finally calls are annotated to include if they are within a defined exon (exons
boundaries were extracted by the R package annotatr, which uses the
`TxDb.Mmusculus.UCSC.mm10.knownGene` resource.

Annotations in that package were drawn from resources at
UCSC on _2019-10-21 20:52:26 +0000 (Mon, 21 Oct 2019)_ and based on the
_mm10_ genome based on the knownGene table)

Structural variant type are classed into the following types:
    - _INS_ – Insertion
    - _INV_ – Inversion
    - _DEL_ – Deletion
    - _DUP_ – Duplication
    - _TRA_ – Translocation
