#!/usr/bin/env python3
"""
Select nascent read names directly from deduplicated BAM files.

This is a storage-efficient replacement for:
  1. bam_to_align.sh, which writes a huge sam2tsv base table
  2. select_nascent_reads.R, which reads that table into R

Usage:
    python3 select_nascent_reads_pysam.py DEDUP_BAM_FOLDER SAMPLE_ID_FILE READ_NAMES_OUTPUT_FOLDER CORE SNP_TABLE REF_FASTA [NASCENT_BAM_OUTPUT_FOLDER]

The CORE argument is accepted for compatibility with the old R script call, but
this implementation streams one BAM at a time in a single process.

If NASCENT_BAM_OUTPUT_FOLDER is provided, the script also writes:
  <sample>_nascent_dedup.bam
  <sample>_nascent_dedup_sorted.bam
  <sample>_nascent_dedup_sorted.bam.bai
"""

import csv
import os
import sys
from typing import Dict, Iterable, Optional, Set, Tuple

import pysam


QUALITY_FILTER_ASCII = 45
TARGET_MUT_FILTER_RATE = 0.3
PROGRESS_EVERY = 1_000_000

CIGAR_M = 0
CIGAR_I = 1
CIGAR_D = 2
CIGAR_N = 3
CIGAR_S = 4
CIGAR_H = 5
CIGAR_P = 6
CIGAR_EQ = 7
CIGAR_X = 8

CONSUMES_QUERY = {CIGAR_M, CIGAR_I, CIGAR_S, CIGAR_EQ, CIGAR_X}
CONSUMES_REF = {CIGAR_M, CIGAR_D, CIGAR_N, CIGAR_EQ, CIGAR_X}
COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def read_sample_names(path: str) -> list:
    names = []
    with open(path, "r", newline="") as handle:
        for row in csv.reader(handle):
            for value in row:
                value = value.strip()
                if value:
                    names.append(value)
    return names


def read_snp_keys(path: str) -> Set[Tuple[str, int, str, str]]:
    snps = set()
    with open(path, "r", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"Chrom", "Position", "Ref", "Var"}
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(f"SNP table is missing required columns: {', '.join(sorted(missing))}")

        for row in reader:
            try:
                pos = int(row["Position"])
            except ValueError:
                continue
            snps.add((
                row["Chrom"],
                pos,
                row["Ref"].upper(),
                row["Var"].upper(),
            ))
    return snps


def plus_strand_query_base(read: pysam.AlignedSegment, query_pos: int) -> str:
    base = read.query_sequence[query_pos]
    if read.is_reverse:
        return base.translate(COMPLEMENT).upper()
    return base.upper()


def is_target_mutation(read: pysam.AlignedSegment, ref_base: str, read_base: str) -> bool:
    if read.is_reverse:
        return ref_base == "A" and read_base == "G"
    return ref_base == "T" and read_base == "C"


def iter_m_bases(read: pysam.AlignedSegment, fasta: pysam.FastaFile) -> Iterable[Tuple[int, str, str, int]]:
    query_pos = 0
    ref_pos = read.reference_start

    for op, length in read.cigartuples or []:
        if op == CIGAR_M:
            contig = read.reference_name
            ref_seq = fasta.fetch(contig, ref_pos, ref_pos + length).upper()
            for offset, ref_base in enumerate(ref_seq):
                qpos = query_pos + offset
                read_base = plus_strand_query_base(read, qpos)
                base_qual = read.query_qualities[qpos] if read.query_qualities is not None else None
                yield ref_pos + offset, ref_base, read_base, base_qual

        if op in CONSUMES_QUERY:
            query_pos += length
        if op in CONSUMES_REF:
            ref_pos += length


def select_sample(
    sample_name: str,
    bam_folder: str,
    read_names_output_folder: str,
    snps: Set[Tuple[str, int, str, str]],
    fasta: pysam.FastaFile,
    nascent_bam_output_folder: Optional[str] = None,
) -> float:
    input_bam = os.path.join(bam_folder, f"{sample_name}_dedup_sorted.bam")
    output_path = os.path.join(read_names_output_folder, f"{sample_name}.txt")

    print(f"Process sample:  {sample_name}", file=sys.stderr, flush=True)

    if not os.path.exists(input_bam) or os.path.getsize(input_bam) == 0:
        open(output_path, "w").close()
        if nascent_bam_output_folder is not None:
            write_empty_bam(sample_name, input_bam, nascent_bam_output_folder)
        return -1.0

    total_reads = 0
    output_reads = 0
    min_phred_quality = QUALITY_FILTER_ASCII - 33

    output_bam_handle = None
    output_bam_path = None
    if nascent_bam_output_folder is not None:
        os.makedirs(nascent_bam_output_folder, exist_ok=True)
        output_bam_path = os.path.join(nascent_bam_output_folder, f"{sample_name}_nascent_dedup.bam")

    with pysam.AlignmentFile(input_bam, "rb") as bam, open(output_path, "w") as out_handle:
        if output_bam_path is not None:
            output_bam_handle = pysam.AlignmentFile(output_bam_path, "wb", template=bam)

        for read in bam.fetch(until_eof=True):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            if read.is_paired:
                continue
            if read.is_reverse:
                flag_matches_r_script = read.flag == 16
            else:
                flag_matches_r_script = read.flag == 0
            if not flag_matches_r_script:
                continue
            if read.query_sequence is None:
                continue

            total_reads += 1
            if total_reads % PROGRESS_EVERY == 0:
                print(
                    f"{sample_name}: scanned {total_reads:,} reads, selected {output_reads:,}",
                    file=sys.stderr,
                    flush=True,
                )

            mut_num = 0
            target_mut_num = 0
            chrom = read.reference_name

            for ref_pos0, ref_base, read_base, base_qual in iter_m_bases(read, fasta):
                if ref_base == "." or read_base == "." or read_base == ref_base:
                    continue
                if base_qual is None or base_qual <= min_phred_quality:
                    continue

                ref_pos1 = ref_pos0 + 1
                if (chrom, ref_pos1, ref_base, read_base) in snps:
                    continue

                mut_num += 1
                if is_target_mutation(read, ref_base, read_base):
                    target_mut_num += 1

            if mut_num > 0 and target_mut_num > 0 and (target_mut_num / mut_num) >= TARGET_MUT_FILTER_RATE:
                out_handle.write(read.query_name + "\n")
                if output_bam_handle is not None:
                    output_bam_handle.write(read)
                output_reads += 1

    if output_bam_handle is not None:
        output_bam_handle.close()
        sort_and_index_bam(sample_name, output_bam_path, nascent_bam_output_folder)

    if total_reads == 0:
        open(output_path, "w").close()
        if nascent_bam_output_folder is not None:
            write_empty_bam(sample_name, input_bam, nascent_bam_output_folder)
        return -1.0

    return output_reads / total_reads


def sort_and_index_bam(sample_name: str, output_bam_path: str, nascent_bam_output_folder: str) -> None:
    sorted_bam_path = os.path.join(nascent_bam_output_folder, f"{sample_name}_nascent_dedup_sorted.bam")
    pysam.sort("-o", sorted_bam_path, output_bam_path)
    pysam.index(sorted_bam_path)


def write_empty_bam(sample_name: str, input_bam: str, nascent_bam_output_folder: str) -> None:
    os.makedirs(nascent_bam_output_folder, exist_ok=True)
    output_bam_path = os.path.join(nascent_bam_output_folder, f"{sample_name}_nascent_dedup.bam")

    if os.path.exists(input_bam) and os.path.getsize(input_bam) > 0:
        with pysam.AlignmentFile(input_bam, "rb") as bam:
            with pysam.AlignmentFile(output_bam_path, "wb", template=bam):
                pass
    else:
        open(output_bam_path, "wb").close()

    if os.path.getsize(output_bam_path) > 0:
        sort_and_index_bam(sample_name, output_bam_path, nascent_bam_output_folder)


def write_outputs(read_names_output_folder: str, summary_rows: list) -> None:
    with open(os.path.join(read_names_output_folder, "summary.csv"), "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["sample_name", "mutation_reads_rate"])
        writer.writerows(summary_rows)

    with open(os.path.join(read_names_output_folder, "sample_id.txt"), "w") as handle:
        for sample_name, mut_rate in summary_rows:
            if mut_rate >= 0:
                handle.write(sample_name + "\n")


def main() -> int:
    if len(sys.argv) not in {7, 8}:
        print(__doc__.strip(), file=sys.stderr)
        return 2

    bam_folder, sample_id_file, read_names_output_folder, _core, snp_table, ref_fasta = sys.argv[1:7]
    nascent_bam_output_folder = sys.argv[7] if len(sys.argv) == 8 else None
    os.makedirs(read_names_output_folder, exist_ok=True)
    if nascent_bam_output_folder is not None:
        os.makedirs(nascent_bam_output_folder, exist_ok=True)

    sample_names = read_sample_names(sample_id_file)
    snps = read_snp_keys(snp_table)

    summary_rows = []
    with pysam.FastaFile(ref_fasta) as fasta:
        for sample_name in sample_names:
            mut_rate = select_sample(
                sample_name,
                bam_folder,
                read_names_output_folder,
                snps,
                fasta,
                nascent_bam_output_folder,
            )
            summary_rows.append((sample_name, mut_rate))

    write_outputs(read_names_output_folder, summary_rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
