#!/usr/bin/env python3
"""Split paired FASTQ files into synchronized gzip chunks."""

from __future__ import annotations

import argparse
import gzip
from pathlib import Path
from typing import Iterator, TextIO, Tuple


FastqRecord = Tuple[str, str, str, str]


def open_text(path: str) -> TextIO:
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt", encoding="utf-8")


def read_records(handle: TextIO) -> Iterator[FastqRecord]:
    while True:
        header = handle.readline()
        if not header:
            return
        seq = handle.readline()
        plus = handle.readline()
        qual = handle.readline()
        if not seq or not plus or not qual:
            raise ValueError("Input FASTQ ended in the middle of a record")
        if not header.startswith("@") or not plus.startswith("+"):
            raise ValueError(f"Malformed FASTQ record near header: {header.strip()}")
        yield header, seq, plus, qual


def read_id(header: str) -> str:
    """Return a normalized read ID without common paired-end suffixes."""
    token = header.strip().split()[0]
    if not token.startswith("@"):
        raise ValueError(f"Malformed FASTQ header: {header.strip()}")
    token = token[1:]
    for suffix in ("/1", "/2"):
        if token.endswith(suffix):
            token = token[: -len(suffix)]
    return token


def write_record(handle: gzip.GzipFile, record: FastqRecord) -> None:
    for field in record:
        handle.write(field.encode("utf-8"))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Split paired-end FASTQ files into synchronized gzip chunks."
    )
    parser.add_argument("--r1", required=True)
    parser.add_argument("--r2", required=True)
    parser.add_argument("--chunk-size", required=True, type=int)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--prefix", default="chunk")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.chunk_size <= 0:
        raise SystemExit("--chunk-size must be greater than 0")

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    manifest_path = outdir / "chunks.tsv"

    chunk_index = 0
    pairs_in_chunk = 0
    total_pairs = 0
    r1_out = None
    r2_out = None
    manifest_rows = []

    def close_chunk() -> None:
        nonlocal r1_out, r2_out
        if r1_out is not None:
            r1_out.close()
        if r2_out is not None:
            r2_out.close()
        r1_out = None
        r2_out = None

    try:
        with open_text(args.r1) as r1_handle, open_text(args.r2) as r2_handle:
            r1_records = read_records(r1_handle)
            r2_records = read_records(r2_handle)

            while True:
                try:
                    r1_record = next(r1_records)
                except StopIteration:
                    try:
                        next(r2_records)
                    except StopIteration:
                        break
                    raise ValueError("R2 contains more records than R1")

                try:
                    r2_record = next(r2_records)
                except StopIteration as exc:
                    raise ValueError("R1 contains more records than R2") from exc

                r1_id = read_id(r1_record[0])
                r2_id = read_id(r2_record[0])
                if r1_id != r2_id:
                    raise ValueError(
                        f"Read ID mismatch at pair {total_pairs + 1}: R1={r1_id} R2={r2_id}"
                    )

                if pairs_in_chunk == 0:
                    chunk_index += 1
                    chunk_id = f"{args.prefix}_{chunk_index:06d}"
                    r1_path = outdir / f"{chunk_id}_R1.fastq.gz"
                    r2_path = outdir / f"{chunk_id}_R2.fastq.gz"
                    r1_out = gzip.open(r1_path, "wb")
                    r2_out = gzip.open(r2_path, "wb")
                    manifest_rows.append([chunk_id, str(r1_path), str(r2_path), 0])

                write_record(r1_out, r1_record)
                write_record(r2_out, r2_record)
                pairs_in_chunk += 1
                total_pairs += 1
                manifest_rows[-1][3] = pairs_in_chunk

                if pairs_in_chunk >= args.chunk_size:
                    close_chunk()
                    pairs_in_chunk = 0
    finally:
        close_chunk()

    if total_pairs == 0:
        raise SystemExit("No read pairs found in input FASTQ files")

    with manifest_path.open("w", encoding="utf-8") as manifest:
        manifest.write("chunk_id\tr1_path\tr2_path\tn_read_pairs\n")
        for row in manifest_rows:
            manifest.write("\t".join(str(value) for value in row) + "\n")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
