# HiC-Nap

HiC-Nap is a conservative, restartable, chunk-based Hi-C preprocessing pipeline for paired-end FASTQ files on a local Linux workstation.

This first version stops at per-chunk selected sorted `pairsam.gz` files. Global merging, deduplication, `.pairs.gz`, `.cool`, `.mcool`, and `.hic` generation are intentionally left for a later phase.

## Requirements

The pipeline is written as Bash plus a small Python standard-library FASTQ splitter. The following commands must be available in `PATH`:

- `fastqc`
- `trim_galore`
- `bwa`
- `samtools`
- `pairtools`
- `cooler`
- `gzip`
- `python3`

`bgzip` is not required in this phase.

## Usage

```bash
./hic_chunk_pipeline.sh \
  --sample SAMPLE_ID \
  --r1 /path/to/sample_R1.fastq.gz \
  --r2 /path/to/sample_R2.fastq.gz \
  --genome-name mm10 \
  --genome-fa /path/to/genome.fa \
  --enzyme MboI \
  --workdir /path/to/workdir \
  --outdir /path/to/outdir
```

Optional arguments:

```bash
--threads 14
--chunk-size 10000000
--min-mapq 30
--short-cis-cutoff 2000
--max-chunks 0
--force-init
--status-only
```

Example for a small test run:

```bash
./hic_chunk_pipeline.sh \
  --sample TEST \
  --r1 test_R1.fastq.gz \
  --r2 test_R2.fastq.gz \
  --genome-name mm10 \
  --genome-fa /path/to/mm10.fa \
  --enzyme MboI \
  --workdir /Storage/hic_work \
  --outdir /data/hic_project \
  --threads 4 \
  --chunk-size 25000
```

## Directory Layout

`outdir` stores durable outputs, logs, reference side files, and status:

```text
outdir/
  fastqc/
  logs/SAMPLE/
  genome/
  status/
    .locks/
    SAMPLE/
  chunks/SAMPLE/
    fastq/
    processed/
```

`workdir` stores heavy temporary chunk intermediates:

```text
workdir/
  SAMPLE/
    tmp/
    chunks/
```

Each completed chunk produces:

```text
outdir/chunks/SAMPLE/processed/chunk_000001/chunk.selected.sorted.pairsam.gz
```

## Resume Behavior

Status files contain exactly one value: `null`, `running`, `done`, or `failed`.

A step is complete only when its status is `done` and its expected output validates. On restart, any missing, invalid, `running`, `failed`, or `null` step is rerun conservatively. The script deletes that step’s output and downstream intermediates before rerunning.

The sample lock is stored at:

```text
outdir/status/.locks/SAMPLE.lock
```

This lock is acquired before `--force-init`, so force reinitialization cannot delete another active sample lock.

Reference preparation is also locked per genome/enzyme under `outdir/genome/.locks/`, so multiple samples can safely share reference preparation.

## Nightly Partial Runs

Use `--max-chunks` to process only a limited number of unfinished chunks in one invocation:

```bash
./hic_chunk_pipeline.sh \
  --sample SAMPLE_ID \
  --r1 /path/to/sample_R1.fastq.gz \
  --r2 /path/to/sample_R2.fastq.gz \
  --genome-name mm10 \
  --genome-fa /path/to/genome.fa \
  --enzyme MboI \
  --workdir /path/to/workdir \
  --outdir /path/to/outdir \
  --max-chunks 2
```

If the script stops cleanly because the chunk limit is reached and chunks remain incomplete, `pipeline.status` is reset to `null`. Errors and interrupts set `pipeline.status` to `failed`.

## Status Summary

```bash
./hic_chunk_pipeline.sh \
  --sample SAMPLE_ID \
  --r1 /path/to/sample_R1.fastq.gz \
  --r2 /path/to/sample_R2.fastq.gz \
  --genome-name mm10 \
  --genome-fa /path/to/genome.fa \
  --enzyme MboI \
  --workdir /path/to/workdir \
  --outdir /path/to/outdir \
  --status-only
```

`--status-only` prints global and per-step completion counts without running processing.

## Current Non-goals

This version does not implement:

- global `pairtools merge`
- global `pairtools dedup`
- `.pairs.gz` output
- pairix indexing
- `cooler cload`
- `cooler zoomify`
- Juicer `.hic` generation
- `.mcool` output
- Hi-C QC summary
- preseq
- BAM output
