#!/usr/bin/env python3
"""
Select nascent read names directly from deduplicated BAM files in parallel.

Usage:
    python3 select_nascent_reads_pysam_parallel.py DEDUP_BAM_FOLDER SAMPLE_ID_FILE READ_NAMES_OUTPUT_FOLDER CORE SNP_TABLE REF_FASTA [NASCENT_BAM_OUTPUT_FOLDER]

This keeps the same inputs and outputs as select_nascent_reads_pysam.py, but
uses CORE worker processes to process samples concurrently. Each worker opens
its own reference FASTA handle because pysam file handles should not be shared
between processes.
"""

import os
import sys
from concurrent.futures import ProcessPoolExecutor
from typing import Optional, Set, Tuple

import pysam

from select_nascent_reads_pysam import read_sample_names, read_snp_keys, select_sample, write_outputs


_BAM_FOLDER: Optional[str] = None
_READ_NAMES_OUTPUT_FOLDER: Optional[str] = None
_SNPS: Optional[Set[Tuple[str, int, str, str]]] = None
_REF_FASTA: Optional[str] = None
_NASCENT_BAM_OUTPUT_FOLDER: Optional[str] = None
_FASTA: Optional[pysam.FastaFile] = None


def parse_core(core: str, sample_count: int) -> int:
    try:
        workers = int(core)
    except ValueError as exc:
        raise ValueError(f"CORE must be an integer, got {core!r}") from exc

    if workers < 1:
        raise ValueError(f"CORE must be >= 1, got {workers}")

    return min(workers, max(sample_count, 1))


def init_worker(
    bam_folder: str,
    read_names_output_folder: str,
    snps: Set[Tuple[str, int, str, str]],
    ref_fasta: str,
    nascent_bam_output_folder: Optional[str],
) -> None:
    global _BAM_FOLDER
    global _READ_NAMES_OUTPUT_FOLDER
    global _SNPS
    global _REF_FASTA
    global _NASCENT_BAM_OUTPUT_FOLDER
    global _FASTA

    _BAM_FOLDER = bam_folder
    _READ_NAMES_OUTPUT_FOLDER = read_names_output_folder
    _SNPS = snps
    _REF_FASTA = ref_fasta
    _NASCENT_BAM_OUTPUT_FOLDER = nascent_bam_output_folder
    _FASTA = pysam.FastaFile(ref_fasta)


def select_sample_worker(sample_name: str) -> Tuple[str, float]:
    if (
        _BAM_FOLDER is None
        or _READ_NAMES_OUTPUT_FOLDER is None
        or _SNPS is None
        or _REF_FASTA is None
        or _FASTA is None
    ):
        raise RuntimeError("Worker was not initialized")

    mut_rate = select_sample(
        sample_name,
        _BAM_FOLDER,
        _READ_NAMES_OUTPUT_FOLDER,
        _SNPS,
        _FASTA,
        _NASCENT_BAM_OUTPUT_FOLDER,
    )
    return sample_name, mut_rate


def main() -> int:
    if len(sys.argv) not in {7, 8}:
        print(__doc__.strip(), file=sys.stderr)
        return 2

    bam_folder, sample_id_file, read_names_output_folder, core, snp_table, ref_fasta = sys.argv[1:7]
    nascent_bam_output_folder = sys.argv[7] if len(sys.argv) == 8 else None

    os.makedirs(read_names_output_folder, exist_ok=True)
    if nascent_bam_output_folder is not None:
        os.makedirs(nascent_bam_output_folder, exist_ok=True)

    sample_names = read_sample_names(sample_id_file)
    snps = read_snp_keys(snp_table)
    workers = parse_core(core, len(sample_names))

    if workers == 1:
        with pysam.FastaFile(ref_fasta) as fasta:
            summary_rows = [
                (
                    sample_name,
                    select_sample(
                        sample_name,
                        bam_folder,
                        read_names_output_folder,
                        snps,
                        fasta,
                        nascent_bam_output_folder,
                    ),
                )
                for sample_name in sample_names
            ]
    else:
        print(f"Processing {len(sample_names)} samples with {workers} workers", file=sys.stderr, flush=True)
        with ProcessPoolExecutor(
            max_workers=workers,
            initializer=init_worker,
            initargs=(bam_folder, read_names_output_folder, snps, ref_fasta, nascent_bam_output_folder),
        ) as executor:
            summary_rows = list(executor.map(select_sample_worker, sample_names))

    write_outputs(read_names_output_folder, summary_rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
