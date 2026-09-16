# Nascent Read Processing

This copy adds an optional nascent GEX branch to the PerturbSci pipeline without changing the original total GEX, GDO, or Seurat steps.

Set `run_nascent=true` in `main_pipeline.sh` to enable it. The required extra inputs are:

- `ref_genome_fa`: reference FASTA used by JVarkit `sam2tsv`.
- `ref_SNP_var_file`: bulk/background SNP table with `Chrom`, `Position`, `Ref`, and `Var` columns.
- `sam2tsv_cmd`: JVarkit `sam2tsv` command, or a path to `sam2tsv.jar`.
- `picard_cmd`: Picard command, or a path to `picard.jar`.
- `java_bin`: Java executable.

The PerturbSci read roles are:

- `R1`: UMI plus RT barcode.
- `R2`: ligation barcode.
- `R3`: transcript read.

`barcode_extraction.py` already follows that layout. It writes the transcript sequence to `BC_attach/<sample>.R2.fastq.gz` and places the corrected cell barcode plus UMI in the read name in the format expected by `umi_tools`.

The nascent branch runs after `umi_tools dedup` creates `dedup/<sample>_dedup_sorted.bam`:

1. `bam_to_align.sh` converts the deduplicated BAM to `nascent_align/<sample>.align` with JVarkit `sam2tsv`.
2. `select_nascent_reads.R` removes known SNPs and selects reads enriched for SLAM-style conversions:
   - plus-strand reads: `T>C`
   - minus-strand reads: `A>G`
   - base quality must be greater than 45
   - target conversion fraction must be at least 0.3
3. `filter_nascent_bam.sh` uses Picard `FilterSamReads` to keep only selected reads from the deduplicated BAM.
4. `nascent_count.sh` runs `umi_tools count` with the same per-cell/per-gene settings as the total GEX matrix.

New outputs are written under `gex_processing`:

- `nascent_align/`: per-sample `sam2tsv` tables.
- `nascent_read_names/`: selected nascent read IDs and `summary.csv`.
- `nascent_dedup/`: nascent-only deduplicated BAMs.
- `nascent_count/`: nascent cell-by-gene count matrices.
