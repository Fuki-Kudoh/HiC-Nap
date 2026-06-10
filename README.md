# HiC-Nap

HiC-Nap is a conservative, restartable, chunk-based Hi-C preprocessing pipeline for paired-end FASTQ files on a local Linux workstation.

Version 0.2.x runs chunk-level preprocessing, then globally merges and deduplicates all selected chunk `pairsam.gz` files to produce a BGZF-compressed, pairix-indexed sample-level valid pairs file:

```text
outdir/pairs/SAMPLE.valid.pairs.gz
outdir/pairs/SAMPLE.valid.pairs.gz.px2
```

Version 0.3.0 adds a restartable matrix-generation phase after valid pairs are complete.

## Requirements

The pipeline is written as Bash plus a small Python standard-library FASTQ splitter. The following commands must be available in `PATH`:

- `fastqc`
- `trim_galore`
- `bwa`
- `samtools`
- `pairtools`
- `pairix`
- `bgzip`
- `cooler`
- `juicer_tools`
- `gzip`
- `python3`

HiC-Nap calls `bgzip` during `pairtools split` to ensure the final valid pairs output is compatible with pairix indexing.

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
--chunk-validate-mode light|strict
--min-mapq 30
--short-cis-cutoff 2000
--max-chunks 0
--stop-after-chunks
--skip-chunks
--merge-memory 4G
--merge-max-nmerge 8
--force-init
--status-only
--no-matrix
--no-cool
--no-hic
--matrix-only
--cool-base-res 10000
--mcool-resolutions 10000,50000,100000,250000,500000,1000000
```

Every `--mcool-resolutions` value must be an integer multiple of `--cool-base-res`. For 25 kb output, use a 5 kb base resolution:

```bash
--cool-base-res 5000
--mcool-resolutions 5000,10000,25000,50000,100000,250000,500000,1000000
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
  pairs/
  stats/
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
    merge/
```

Each completed chunk produces:

```text
outdir/chunks/SAMPLE/processed/chunk_000001/chunk.selected.sorted.pairsam.gz
```

The completed sample-level phase produces BGZF-compressed and pairix-indexed valid pairs:

```text
outdir/pairs/SAMPLE.valid.pairs.gz
outdir/pairs/SAMPLE.valid.pairs.gz.px2
outdir/stats/SAMPLE.dedup.stats
outdir/status/SAMPLE/chunk_pairsam.list
```

The v0.3.0 matrix phase then produces:

```text
outdir/cool/SAMPLE.10000.cool
outdir/cool/SAMPLE.mcool
outdir/hic/SAMPLE.hic
```

## Resume Behavior

Status files contain exactly one value: `null`, `running`, `done`, or `failed`.

A step is complete only when its status is `done` and its expected output validates. On restart, any missing, invalid, `running`, `failed`, or `null` step is rerun conservatively. The script deletes that step’s output and downstream intermediates before rerunning.

In v0.2.0, `pipeline.status = done` means chunk outputs validate and the final valid pairs file plus pairix index exist and validate. Completed chunks alone are not enough for final pipeline completion.

In v0.3.0 default runs, `pipeline.status = done` means valid pairs and enabled matrix outputs validate. Use `--no-matrix` to stop after the v0.2.x valid-pairs phase.

The sample lock is stored at:

```text
outdir/status/.locks/SAMPLE.lock
```

This lock is acquired before `--force-init`, so force reinitialization cannot delete another active sample lock.

Reference preparation is also locked per genome/enzyme under `outdir/genome/.locks/`, so multiple samples can safely share reference preparation.

## Chunk Manifest Validation

HiC-Nap validates existing FASTQ chunks before resuming from a previous run.

By default, validation uses `--chunk-validate-mode light`, which checks that each chunk has a non-empty ID, is listed with a positive `n_read_pairs` value in `chunks.tsv`, and has non-empty R1/R2 FASTQs that pass `gzip -t`. This avoids re-counting every read on restart.

For maximum conservativeness, use:

```bash
--chunk-validate-mode strict
```

Strict mode decompresses all chunk FASTQs and verifies that R1/R2 read counts match the `n_read_pairs` column in `chunks.tsv`. This is safer but can be slow for large Hi-C datasets.

In light mode, a gzip-valid FASTQ with the wrong read count may not be detected during manifest validation. It should still fail later when that chunk is processed by `trim_galore` or downstream validation.

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

## Chunk-only Runs

Use `--stop-after-chunks` to reproduce the v0.1.x behavior and stop after chunk preprocessing:

```bash
./hic_chunk_pipeline.sh \
  --sample SAMPLE_ID \
  --r1 /path/to/sample_R1.fastq.gz \
  --r2 /path/to/sample_R2.fastq.gz \
  --genome-name mm10 \
  --genome-fa /path/to/mm10.fa \
  --enzyme MboI \
  --workdir /path/to/workdir \
  --outdir /path/to/outdir \
  --stop-after-chunks
```

This leaves `pipeline.status` as `null` and does not create `SAMPLE.valid.pairs.gz`.

## Continuing from v0.1.x chunk outputs

If you already processed a sample with HiC-Nap v0.1.x and have completed chunk-level outputs, you can continue directly to sample-level pairs generation:

```bash
./hic_chunk_pipeline.sh \
  --sample SAMPLE_ID \
  --r1 /path/to/sample_R1.fastq.gz \
  --r2 /path/to/sample_R2.fastq.gz \
  --genome-name mm10 \
  --genome-fa /path/to/mm10.fa \
  --enzyme MboI \
  --workdir /path/to/workdir \
  --outdir /path/to/outdir \
  --skip-chunks
```

This validates the existing `chunk.selected.sorted.pairsam.gz` files and then runs global merge, global deduplication, pairs output, and pairix indexing.

## Matrix Generation

By default, HiC-Nap v0.3.0 runs matrix generation after `outdir/pairs/SAMPLE.valid.pairs.gz` and its `.px2` pairix index validate. If those valid-pairs outputs already exist and validate, the script skips chunk preprocessing, merge, deduplication, split-pairs, and pairix indexing, then proceeds directly to matrix generation.

Matrix generation writes restartable statuses:

```text
outdir/status/SAMPLE/cool.status
outdir/status/SAMPLE/mcool.status
outdir/status/SAMPLE/hic.status
outdir/status/SAMPLE/matrix.status
```

The matrix commands are:

```bash
cooler cload pairix -p THREADS outdir/genome/GENOME.chrom.sizes:COOL_BASE_RES outdir/pairs/SAMPLE.valid.pairs.gz outdir/cool/SAMPLE.COOL_BASE_RES.cool
cooler zoomify -p THREADS --balance -r MCOOL_RESOLUTIONS -o outdir/cool/SAMPLE.mcool outdir/cool/SAMPLE.COOL_BASE_RES.cool
juicer_tools pre -j THREADS outdir/pairs/SAMPLE.valid.pairs.gz outdir/hic/SAMPLE.hic outdir/genome/GENOME.chrom.sizes
```

The `.hic` command uses the `chrom.sizes` path rather than only a genome name, so custom genomes are supported.

To upgrade an existing v0.2.x output directory directly to v0.3.0 matrices:

```bash
./hic_chunk_pipeline.sh \
  --sample SAMPLE \
  --genome-name GENOME \
  --genome-fa GENOME.fa \
  --outdir OUTDIR \
  --matrix-only
```

Matrix-only mode requires existing valid pairs and pairix index outputs, creates or validates `outdir/genome/GENOME.chrom.sizes` from `GENOME.fa.fai`, and does not require R1, R2, enzyme, workdir, a restriction digest, or a BWA index.

## Progress Output

The pipeline prints line-based progress messages to stderr and appends them to `outdir/logs/SAMPLE/global.log`, so long runs can be followed from a terminal, tmux pane, SSH session, or redirected job log.

Example output:

```text
[HiC-Nap] run summary:
[HiC-Nap]   sample: SAMPLE_ID
[HiC-Nap]   input R1: /path/to/sample_R1.fastq.gz
[HiC-Nap]   input R2: /path/to/sample_R2.fastq.gz
[HiC-Nap]   genome name: mm10
[HiC-Nap]   enzyme: MboI
[HiC-Nap]   threads: 14
[HiC-Nap]   chunk size: 10000000
[HiC-Nap]   chunk validate mode: light
[HiC-Nap]   max chunks: 0
[HiC-Nap]   stop after chunks: 0
[HiC-Nap]   skip chunks: 0
[HiC-Nap]   merge memory: 4G
[HiC-Nap]   merge max nmerge: 8
[HiC-Nap]   workdir: /path/to/workdir
[HiC-Nap]   outdir: /path/to/outdir
[HiC-Nap] chunk splitting complete:
[HiC-Nap]   total chunks: 42
[HiC-Nap]   total read pairs: 420000000
[HiC-Nap]   chunk size: 10000000
[HiC-Nap] chunk 3/42 chunk_000003: trim_galore
[HiC-Nap] chunk 3/42 chunk_000003: bwa_mem
[HiC-Nap] chunk 3/42 chunk_000003: parse
[HiC-Nap] chunk 3/42 chunk_000003: sort
[HiC-Nap] chunk 3/42 chunk_000003: restrict
[HiC-Nap] chunk 3/42 chunk_000003: select
[HiC-Nap] chunk 3/42 chunk_000003: done
[HiC-Nap] chunk 7/42 chunk_000007: already complete, skipping
[HiC-Nap] starting sample-level pairs generation
[HiC-Nap] validating chunk-level selected pairsam outputs
[HiC-Nap] validated chunk outputs: 42/42
[HiC-Nap] merging chunk pairsam files
[HiC-Nap] global deduplication
[HiC-Nap] writing valid pairs file
[HiC-Nap] indexing valid pairs with pairix
[HiC-Nap] valid pairs complete: /path/to/outdir/pairs/SAMPLE_ID.valid.pairs.gz
[HiC-Nap] final summary:
[HiC-Nap]   total chunks: 42
[HiC-Nap]   completed chunks: 42
[HiC-Nap]   incomplete chunks: 0
[HiC-Nap]   final pipeline.status: done
[HiC-Nap]   valid pairs.status: done
[HiC-Nap]   status directory: /path/to/outdir/status/SAMPLE_ID
[HiC-Nap]   logs directory: /path/to/outdir/logs/SAMPLE_ID
```

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

The status summary includes sample-level pairs statuses:

```text
Sample-level pairs:
  chunk_outputs: null/running/done/failed
  merge: null/running/done/failed
  dedup: null/running/done/failed
  split_pairs: null/running/done/failed
  pairix: null/running/done/failed
  valid_pairs: null/running/done/failed

Outputs:
  valid pairs: outdir/pairs/SAMPLE.valid.pairs.gz
  pairix index: outdir/pairs/SAMPLE.valid.pairs.gz.px2
```

## Current Non-goals

This version does not implement:

- `cooler cload`
- `cooler zoomify`
- `.cool` output
- `.mcool` output
- Juicer `.hic` generation
- Hi-C QC summary
- preseq
- BAM output
- multi-sample workflow
- parallel sample scheduling
