<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/mochi_logo_dark.svg">
    <img src="docs/mochi_logo_light.svg" width="420" alt="MOCHI">
  </picture>
</p>

<p align="center"><b>M</b>apping <b>O</b>f <b>C</b>ryptic <b>H</b>idden <b>I</b>ntronic transcripts</p>

MOCHI is a Nextflow pipeline for detecting candidate intronic RNA units from stranded short-read RNA-seq data, combining reference-free k-mer analysis with alignment-based coverage evidence.

<p align="center">
  <img src="docs/WorkflowDiagram.svg" width="800" alt="MOCHI workflow">
</p>

### How it works

1. **k-mer counting.** Reads are trimmed and k-mers are counted per sample with [Jellyfish](https://github.com/gmarcais/Jellyfish), then joined into a k-mer count matrix across samples.
2. **Masking and assembly.** k-mers present in the annotated transcriptome are removed with [KaMRaT](https://github.com/Transipedia/KaMRaT); the remaining k-mers are merged into contigs.
3. **Contig annotation.** Contigs are aligned to the genome and to a repeat reference with [Annotate-contigs](https://github.com/Transipedia/Annotate-contigs); contigs located in introns are kept.
4. **Alignment-based filtering.** Reads are aligned with [STAR](https://github.com/alexdobin/STAR). Introns that are spliced, retained, or overlap exon extensions are excluded.
5. **Candidate reporting.** The remaining introns are summarised per intron, optionally tested for differential expression (DESeq2), and filtered by read coverage.




## Requirements

- Linux (x86_64)
- [Apptainer](https://apptainer.org/docs/admin/main/installation.html) ≥ 1.1
- ~40 GB RAM (for the alignment step)  other steps are light

Nothing else needs to be installed; Nextflow and all other software are bundled in the container.

## Installation

### Option 1: Download prebuilt container (recommended)

```bash
wget -O mochi.sif "https://zenodo.org/records/23592606/files/mochi.sif?download=1"
md5sum mochi.sif   # compare with the MD5 shown on Zenodo
```

`mochi.sif` contains the complete pipeline and all its dependencies.

### Option 2: Build from source

```bash
git clone https://github.com/dcansugurer/MOCHI.git
cd MOCHI
apptainer build mochi.sif container/mochi.def
apptainer test mochi.sif
```

---

## Usage

A run needs `mochi.sif`, the reference files, paired-end FASTQ files and a metafile. Templates of the run script and metafile can be downloaded and edited:

```bash
wget https://raw.githubusercontent.com/dcansugurer/MOCHI/v1.0.0/examples/run_mochi.sh
wget https://raw.githubusercontent.com/dcansugurer/MOCHI/v1.0.0/examples/metafile_single.txt   # or metafile_DE.txt
```

### Required inputs

| Input | Parameter | Details |
|---|---|---|
| Metafile | `--metafile` | [Metafile](#metafile) |
| Paired FASTQ files | `--raw_reads` | quoted glob, e.g. `'/data/fastq/*{_R1,_R2}.fastq.gz'`; see [FASTQ files](#fastq-files) |
| Genome FASTA | `--genome` | may be gzipped |
| Transcriptome FASTA | `--transcriptome` | **uncompressed** |
| Annotation GTF | `--annotation` | may be gzipped |
| Exon BED | `--exon_file` | [Exon and intron BED files](#exon-and-intron-bed-files) |
| Intron BED | `--intron_file` | [Exon and intron BED files](#exon-and-intron-bed-files) |
| Repeat FASTA | `--repeat` | [Repeat reference](#repeat-reference) |
| Library type | `--lib_type` | `rf` (dUTP; R2 = sense) or `fr` |

### Mode parameters (choose one)

| Mode | Parameters | Metafile |
|---|---|---|
| Candidate discovery | `--sample_name <condition>` | one condition |
| Candidate discovery + Differential expression (DESeq2) | `--DE --control_group <condition> --test_group <condition>` | two conditions, ≥ 2 replicates each |

### Index parameters (choose one)

| Run | Parameters |
|---|---|
| First run | `--index` (builds STAR and minimap2 indices) |
| Later runs | `--human_index <dir> --star_index <dir>/star` (see [Alignment indices](#alignment-indices)) |

### Example

```bash
apptainer run --cleanenv -B /path/to/data mochi.sif \
    --metafile      metafile.txt \
    --raw_reads     '/path/to/data/fastq/*{_R1,_R2}.fastq.gz' \
    --genome        /path/to/data/references/GRCh38.p14.genome.fa.gz \
    --transcriptome /path/to/data/references/gencode.v47.transcripts.fa \
    --annotation    /path/to/data/references/gencode.v47.annotation.gtf.gz \
    --exon_file     /path/to/data/references/gencode_v47_exons.bed \
    --intron_file   /path/to/data/references/gencode_v47_introns.bed \
    --repeat        /path/to/data/references/human_repeat_ref.fasta \
    --lib_type      rf \
    --index \
    --sample_name   TEST
```

- `-B /path/to/data` makes the input directory visible inside the container; every input must be under a bound directory (`-B /data,/scratch` for several).
- `--cleanenv` keeps the host environment (e.g. an active conda environment) out of the container.
- Add `-resume` to continue an interrupted run from the last completed step.
- `apptainer run-help mochi.sif` prints a summary of all options.

---

## Input details

### Reference files

>[!IMPORTANT]
> - All reference files must come from the **same genome assembly and annotation release**.
> - The transcriptome FASTA must be **uncompressed**. A gzipped file is not rejected but silently masks only a small fraction of annotated k-mers, which changes the candidate set.
> - BED files must be uncompressed, sorted with `sort -k1,1V -k2,2n -k3,3n`, contain strand information (6 columns), and include **canonical chromosomes only**. Scaffolds and alternative contigs cause the pipeline to fail.

### Exon and intron BED files

**GENCODE v47 / GRCh38:** ready-to-use files are available in the [Zenodo record](https://doi.org/10.5281/zenodo.23592606):

```bash
wget -O - "https://zenodo.org/records/23592606/files/sorted_gencode_v47_introns_filtered.bed?download=1" > gencode_v47_introns.bed
wget -O - "https://zenodo.org/records/23592606/files/sorted_gencode_v47_exons_filtered.bed?download=1" > gencode_v47_exons.bed
```

**Other annotation releases:** Exon and intron BED tables can be retrieved from [UCSC Table Browser](https://genome.ucsc.edu/cgi-bin/hgTables). Non-canonical chromosomes should be removed and BED files should be sorted.

The name column must keep the UCSC format (e.g. `ENST00000456328.2_intron_(flank)0_0_chr1_12228_f`); MOCHI parses transcript IDs and intron numbers from it.

### Repeat reference

Human repeat reference (Dfam), from the [Zenodo record](https://doi.org/10.5281/zenodo.23592606):

```bash
wget -O human_repeat_ref.fasta "https://zenodo.org/records/23592606/files/human_repeat_ref.fasta?download=1"
```

### STAR and minimap2 index files

With `--index`, STAR and minimap2 index files are built on the first run and saved to:

```
output/3_contig_annotation/human_index/
├── star/
└── minimap2/human.mmi
```

On later runs, replace `--index` with `--human_index <path>/human_index --star_index <path>/human_index/star`. Copy `human_index/` to a permanent location before deleting `output/`. Reused indices must match the genome and annotation of the run.

### FASTQ files

Paired-end, **stranded** reads (`rf` = dUTP, or `fr`); unstranded libraries are not supported. In the `--raw_reads` glob, the part matched by `*` is the sample ID and must match the `sample_id` column of the metafile.

### Metafile

Tab-separated, with a header line:

```
condition	sample_id	replicate
CTRL	sample_A	1
CTRL	sample_B	2
TEST	sample_C	1
TEST	sample_D	2
```

Use one condition for candidate discovery (`--sample_name`) and two conditions with at least 2 replicates each for DE (`--control_group`, `--test_group`). Templates: [`metafile_single.txt`](examples/metafile_single.txt), [`metafile_DE.txt`](examples/metafile_DE.txt).

---

## Output

All results are written to `output/` in the directory where MOCHI is launched. Intermediate files are kept so that a run with modified settings can reuse them instead of starting from scratch.

```
output/
├── 1_preprocessed_reads/
│   ├── 1_quality_check/      # FastQC reports
│   ├── 2_trimmed_reads/      # Trim Galore output
│   └── 3_revComp/            # reverse-complemented reads used for k-mer counting
├── 2_kmer_processing/
│   ├── 1_jellyfish/          # per-sample k-mer counts
│   └── 2_kamrat/             # joined, masked and merged k-mer matrices (contigs)
├── 3_contig_annotation/
│   ├── human_index/          # STAR and minimap2 indices (with --index; reusable)
│   └── output/               # Annotate-contigs results
├── 4_aligned_sequences/
│   ├── bams/                 # STAR alignments (sorted BAM + index)
│   └── logs/                 # STAR logs and gene counts
└── 5_candidates/             # final results
```

### Final results (`5_candidates/`)

`<label>` is `<test>_vs_<control>` in DE mode and the `--sample_name` value in single-condition mode.

| File | Description |
|---|---|
| `<label>_covered_candidates.tsv` | **Final candidate intronic RNA units**, after all filters including read coverage |
| `<label>_candidates.tsv` | Candidates before the coverage filter |
| `not_covered_introns.txt` | Candidate introns removed by the coverage filter |
| `retained_introns.bed` | Introns classified as retained (excluded) |
| `intron_coverage_full.tsv` | Strand-specific coverage breadth for all introns |
| `extended_exons.bed` | Exon-boundary extensions supported by reads |
| `extension_containing_introns.txt` | Introns containing exon extensions (excluded) |
| `final_merged_splices.txt` | Splice junctions merged across samples |
| `conditions.Rdata` | Run labels used internally |

### Candidate table columns

| Column | Description |
|---|---|
| `intron_id` | Intron identifier (from the intron BED) |
| `chr`, `start`, `end` | Candidate coordinates. Single-condition mode: span of the intron's contigs. DE mode: contig with the smallest adjusted p-value |
| `intron_number`, `strand`, `transcript_id` | Host intron position, strand and transcript |
| `intron_gene_id`, `intron_gene_name`, `intron_biotype` | Host gene and transcript biotype (from the GTF) |
| `number_of_contigs`, `average_contig_size` | Number and mean length of contigs assigned to the intron |
| `min_contig_pvalue` | DE mode: smallest adjusted p-value (FDR) among the intron's contigs |
| `mean_FC` | DE mode: mean log2 fold change (test vs control) across the intron's contigs |
| `<sample>_mean` | Mean contig count per sample (size-factor normalised in single-condition mode) |
| `is_spliced`, `is_canonical` | Intron overlaps a detected splice junction; host transcript is Ensembl canonical |


## License

TODO: licences
