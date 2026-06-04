#!/usr/bin/env bash
set -Eeuo pipefail

SAMPLE=""
R1=""
R2=""
GENOME_NAME=""
GENOME_FA=""
ENZYME=""
WORKDIR=""
OUTDIR=""
THREADS=14
CHUNK_SIZE=10000000
MIN_MAPQ=30
SHORT_CIS_CUTOFF=2000
MAX_CHUNKS=0
FORCE_INIT=0
STATUS_ONLY=0

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SPLIT_HELPER="${SCRIPT_DIR}/split_pe_fastq.py"
CURRENT_STATUS_FILE=""
LOCK_FD=9
LOCK_DIR=""
LOCK_HELD=0
REF_LOCK_FD=8
REF_LOCK_DIR=""
REF_LOCK_HELD=0
CURRENT_CHUNK_INDEX=0
CURRENT_CHUNK_TOTAL=0

usage() {
  cat <<'EOF'
Usage:
  ./hic_chunk_pipeline.sh \
    --sample SAMPLE_ID \
    --r1 /path/to/sample_R1.fastq.gz \
    --r2 /path/to/sample_R2.fastq.gz \
    --genome-name mm10 \
    --genome-fa /path/to/genome.fa \
    --enzyme MboI \
    --workdir /path/to/workdir \
    --outdir /path/to/outdir

Optional:
  --threads 14
  --chunk-size 10000000
  --min-mapq 30
  --short-cis-cutoff 2000
  --max-chunks 0
  --force-init
  --status-only
EOF
}

log_msg() {
  local message="$1"
  local ts
  ts="$(date '+%Y-%m-%d %H:%M:%S')"
  if [[ -n "${GLOBAL_LOG:-}" ]]; then
    printf '[%s] %s\n' "$ts" "$message" | tee -a "$GLOBAL_LOG" >&2
  else
    printf '[%s] %s\n' "$ts" "$message" >&2
  fi
}

progress_msg() {
  local message="[HiC-Nap] $1"
  if [[ -n "${GLOBAL_LOG:-}" ]]; then
    printf '%s\n' "$message" | tee -a "$GLOBAL_LOG" >&2
  else
    printf '%s\n' "$message" >&2
  fi
}

chunk_progress_msg() {
  local chunk_id="$1"
  local step_message="$2"
  progress_msg "chunk ${CURRENT_CHUNK_INDEX}/${CURRENT_CHUNK_TOTAL} ${chunk_id}: ${step_message}"
}

die() {
  log_msg "ERROR: $*"
  if [[ "$LOCK_HELD" -eq 1 && -n "${PIPELINE_STATUS:-}" && -f "${PIPELINE_STATUS}" ]]; then
    atomic_set_status "$PIPELINE_STATUS" "failed" || true
  fi
  print_final_summary || true
  exit 1
}

on_error() {
  local exit_code=$?
  local line_no=${BASH_LINENO[0]:-unknown}
  if [[ -n "${CURRENT_STATUS_FILE}" && -f "${CURRENT_STATUS_FILE}" ]]; then
    atomic_set_status "$CURRENT_STATUS_FILE" "failed" || true
  fi
  if [[ "$LOCK_HELD" -eq 1 && -n "${PIPELINE_STATUS:-}" && -f "${PIPELINE_STATUS}" ]]; then
    atomic_set_status "$PIPELINE_STATUS" "failed" || true
  fi
  log_msg "Pipeline failed at line ${line_no} with exit code ${exit_code}"
  print_final_summary || true
  exit "$exit_code"
}

on_interrupt() {
  if [[ -n "${CURRENT_STATUS_FILE}" && -f "${CURRENT_STATUS_FILE}" ]]; then
    atomic_set_status "$CURRENT_STATUS_FILE" "failed" || true
  fi
  if [[ "$LOCK_HELD" -eq 1 && -n "${PIPELINE_STATUS:-}" && -f "${PIPELINE_STATUS}" ]]; then
    atomic_set_status "$PIPELINE_STATUS" "failed" || true
  fi
  log_msg "Interrupted; current step was not marked done"
  print_final_summary || true
  exit 130
}

require_command() {
  local cmd="$1"
  command -v "$cmd" >/dev/null 2>&1 || die "Required command not found in PATH: ${cmd}"
}

atomic_set_status() {
  local status_file="$1"
  local status="$2"
  case "$status" in
    null|running|done|failed) ;;
    *) die "Invalid status '${status}' for ${status_file}" ;;
  esac
  mkdir -p "$(dirname "$status_file")"
  printf '%s\n' "$status" > "${status_file}.tmp"
  mv "${status_file}.tmp" "$status_file"
}

get_status() {
  local status_file="$1"
  if [[ ! -f "$status_file" ]]; then
    printf 'null\n'
    return
  fi
  local status
  status="$(tr -d '[:space:]' < "$status_file")"
  case "$status" in
    null|running|done|failed) printf '%s\n' "$status" ;;
    *) printf 'null\n' ;;
  esac
}

init_one_status() {
  local status_file="$1"
  if [[ ! -f "$status_file" ]]; then
    atomic_set_status "$status_file" "null"
  fi
}

init_status() {
  mkdir -p "$STATUS_DIR" "$CHUNK_STATUS_ROOT" "$LOG_DIR" "$CHUNK_LOG_DIR"
  if [[ "$FORCE_INIT" -eq 1 && -d "$STATUS_DIR" ]]; then
    log_msg "CONFIRMED --force-init: deleting existing state directory for sample '${SAMPLE}'"
    rm -rf "$STATUS_DIR"
    mkdir -p "$STATUS_DIR" "$CHUNK_STATUS_ROOT"
  fi
  init_one_status "$PIPELINE_STATUS"
  init_one_status "$QC_STATUS"
  init_one_status "$CHUNK_SPLIT_STATUS"
  init_one_status "$ALL_CHUNKS_STATUS"
  init_one_status "$FINAL_STATUS"
}

validate_gzip() {
  local path="$1"
  [[ -s "$path" ]] || return 1
  gzip -t "$path" >/dev/null 2>&1
}

count_fastq_reads() {
  local path="$1"
  local lines
  if [[ "$path" == *.gz ]]; then
    lines="$(gzip -cd "$path" | awk 'END { print NR }')"
  else
    lines="$(awk 'END { print NR }' "$path")"
  fi
  [[ $((lines % 4)) -eq 0 ]] || return 1
  printf '%s\n' $((lines / 4))
}

validate_fastq_pair() {
  local r1_path="$1"
  local r2_path="$2"
  [[ -s "$r1_path" && -s "$r2_path" ]] || return 1
  if [[ "$r1_path" == *.gz ]]; then validate_gzip "$r1_path" || return 1; fi
  if [[ "$r2_path" == *.gz ]]; then validate_gzip "$r2_path" || return 1; fi
  local r1_count r2_count
  r1_count="$(count_fastq_reads "$r1_path")" || return 1
  r2_count="$(count_fastq_reads "$r2_path")" || return 1
  [[ "$r1_count" == "$r2_count" ]]
}

validate_pairsam_gz() {
  local path="$1"
  [[ -s "$path" ]] || return 1
  validate_gzip "$path" || return 1
  gzip -cd "$path" | awk '
    BEGIN { pairs_format = 0; columns = 0 }
    /^## pairs format/ { pairs_format = 1 }
    /^#columns:/ {
      if ($0 ~ /readID/ && $0 ~ /chrom1/ && $0 ~ /pos1/ &&
          $0 ~ /chrom2/ && $0 ~ /pos2/ && $0 ~ /strand1/ &&
          $0 ~ /strand2/ && $0 ~ /pair_type/) {
        columns = 1
      }
    }
    END { exit ! (pairs_format && columns) }
  '
}

validate_sam() {
  local path="$1"
  [[ -s "$path" ]] || return 1
  grep -q '^@' "$path"
}

fastqc_prefix() {
  local path="$1"
  local base
  base="$(basename "$path")"
  base="${base%.gz}"
  base="${base%.fastq}"
  base="${base%.fq}"
  printf '%s\n' "$base"
}

prepare_dirs() {
  FASTQC_DIR="${OUTDIR}/fastqc"
  LOG_DIR="${OUTDIR}/logs/${SAMPLE}"
  CHUNK_LOG_DIR="${LOG_DIR}/chunks"
  GENOME_DIR="${OUTDIR}/genome"
  STATUS_ROOT="${OUTDIR}/status"
  LOCK_ROOT="${STATUS_ROOT}/.locks"
  STATUS_DIR="${OUTDIR}/status/${SAMPLE}"
  CHUNK_STATUS_ROOT="${STATUS_DIR}/chunks"
  CHUNK_FASTQ_DIR="${OUTDIR}/chunks/${SAMPLE}/fastq"
  CHUNK_PROCESSED_ROOT="${OUTDIR}/chunks/${SAMPLE}/processed"
  SAMPLE_WORKDIR="${WORKDIR}/${SAMPLE}"
  SAMPLE_TMPDIR="${SAMPLE_WORKDIR}/tmp"
  SAMPLE_WORK_CHUNKS="${SAMPLE_WORKDIR}/chunks"
  GLOBAL_LOG="${LOG_DIR}/global.log"
  PIPELINE_STATUS="${STATUS_DIR}/pipeline.status"
  QC_STATUS="${STATUS_DIR}/qc.status"
  CHUNK_SPLIT_STATUS="${STATUS_DIR}/chunk_split.status"
  ALL_CHUNKS_STATUS="${STATUS_DIR}/all_chunks.status"
  FINAL_STATUS="${STATUS_DIR}/final.status"
  CHUNKS_TSV="${STATUS_DIR}/chunks.tsv"
  mkdir -p "$FASTQC_DIR" "$LOG_DIR" "$CHUNK_LOG_DIR" "$GENOME_DIR" "$LOCK_ROOT" "$STATUS_DIR" \
    "$CHUNK_STATUS_ROOT" "$CHUNK_FASTQ_DIR" "$CHUNK_PROCESSED_ROOT" \
    "$SAMPLE_TMPDIR" "$SAMPLE_WORK_CHUNKS"
}

prepare_reference() {
  acquire_reference_lock
  log_msg "Preparing reference files if missing"
  local chrom_sizes="${GENOME_DIR}/${GENOME_NAME}.chrom.sizes"
  local frags_bed="${GENOME_DIR}/${GENOME_NAME}.${ENZYME}.frags.bed"
  if [[ ! -s "${GENOME_FA}.fai" ]]; then
    log_msg "Running samtools faidx for ${GENOME_FA}"
    samtools faidx "$GENOME_FA" >> "$GLOBAL_LOG" 2>&1
  fi
  if [[ ! -s "$chrom_sizes" ]]; then
    cut -f1,2 "${GENOME_FA}.fai" > "${chrom_sizes}.tmp"
    mv "${chrom_sizes}.tmp" "$chrom_sizes"
  fi
  if [[ ! -s "$frags_bed" ]]; then
    cooler digest "$chrom_sizes" "$GENOME_FA" "$ENZYME" > "${frags_bed}.tmp"
    mv "${frags_bed}.tmp" "$frags_bed"
  fi
  local missing_index=0
  local ext
  for ext in amb ann bwt pac sa; do
    if [[ ! -s "${GENOME_FA}.${ext}" ]]; then
      missing_index=1
    fi
  done
  if [[ "$missing_index" -eq 1 ]]; then
    log_msg "BWA index missing; running bwa index for ${GENOME_FA}"
    bwa index "$GENOME_FA" >> "$GLOBAL_LOG" 2>&1
  fi
  release_reference_lock
}

write_step_log_header() {
  local log_file="$1"
  local command_text="$2"
  local output_path="$3"
  {
    printf 'start_time: %s\n' "$(date '+%Y-%m-%d %H:%M:%S')"
    printf 'command: %s\n' "$command_text"
    printf 'output: %s\n' "$output_path"
  } >> "$log_file"
}

write_step_log_footer() {
  local log_file="$1"
  local exit_code="$2"
  {
    printf 'exit_code: %s\n' "$exit_code"
    printf 'end_time: %s\n\n' "$(date '+%Y-%m-%d %H:%M:%S')"
  } >> "$log_file"
}

run_fastqc() {
  local r1_prefix r2_prefix r1_html r1_zip r2_html r2_zip
  r1_prefix="$(fastqc_prefix "$R1")"
  r2_prefix="$(fastqc_prefix "$R2")"
  r1_html="${FASTQC_DIR}/${r1_prefix}_fastqc.html"
  r1_zip="${FASTQC_DIR}/${r1_prefix}_fastqc.zip"
  r2_html="${FASTQC_DIR}/${r2_prefix}_fastqc.html"
  r2_zip="${FASTQC_DIR}/${r2_prefix}_fastqc.zip"
  if [[ "$(get_status "$QC_STATUS")" == "done" && -s "$r1_html" && -s "$r1_zip" && -s "$r2_html" && -s "$r2_zip" ]]; then
    log_msg "FASTQ QC already complete"
    return
  fi
  rm -f "$r1_html" "$r1_zip" "$r2_html" "$r2_zip"
  atomic_set_status "$QC_STATUS" "running"
  CURRENT_STATUS_FILE="$QC_STATUS"
  local log_file="${LOG_DIR}/qc.log"
  local cmd_text="fastqc -t ${THREADS} -o ${FASTQC_DIR} ${R1} ${R2}"
  write_step_log_header "$log_file" "$cmd_text" "$r1_html $r1_zip $r2_html $r2_zip"
  set +e
  fastqc -t "$THREADS" -o "$FASTQC_DIR" "$R1" "$R2" >> "$log_file" 2>&1
  local exit_code=$?
  set -e
  write_step_log_footer "$log_file" "$exit_code"
  CURRENT_STATUS_FILE=""
  if [[ "$exit_code" -eq 0 && -s "$r1_html" && -s "$r1_zip" && -s "$r2_html" && -s "$r2_zip" ]]; then
    atomic_set_status "$QC_STATUS" "done"
  else
    atomic_set_status "$QC_STATUS" "failed"
    die "FASTQ QC failed"
  fi
}

init_chunk_statuses() {
  local chunk_id r1_path r2_path n_read_pairs
  tail -n +2 "$CHUNKS_TSV" | while IFS=$'\t' read -r chunk_id r1_path r2_path n_read_pairs; do
    [[ -n "$chunk_id" ]] || continue
    local chunk_status_dir="${CHUNK_STATUS_ROOT}/${chunk_id}"
    mkdir -p "$chunk_status_dir"
    init_one_status "${chunk_status_dir}/chunk.status"
    init_one_status "${chunk_status_dir}/trim_galore.status"
    init_one_status "${chunk_status_dir}/bwa_mem.status"
    init_one_status "${chunk_status_dir}/parse.status"
    init_one_status "${chunk_status_dir}/sort.status"
    init_one_status "${chunk_status_dir}/restrict.status"
    init_one_status "${chunk_status_dir}/select.status"
  done
}

validate_chunks_manifest() {
  [[ -s "$CHUNKS_TSV" ]] || return 1
  local count
  count="$(tail -n +2 "$CHUNKS_TSV" | awk 'END { print NR }')"
  [[ "$count" -gt 0 ]] || return 1
  local chunk_id r1_path r2_path n_read_pairs observed_pairs
  while IFS=$'\t' read -r chunk_id r1_path r2_path n_read_pairs; do
    [[ "$chunk_id" == "chunk_id" ]] && continue
    [[ -n "$chunk_id" && -s "$r1_path" && -s "$r2_path" ]] || return 1
    validate_fastq_pair "$r1_path" "$r2_path" || return 1
    observed_pairs="$(count_fastq_reads "$r1_path")" || return 1
    [[ "$observed_pairs" == "$n_read_pairs" ]] || return 1
  done < "$CHUNKS_TSV"
}

chunk_total_count() {
  if [[ ! -s "${CHUNKS_TSV:-}" ]]; then
    printf '0\n'
    return
  fi
  tail -n +2 "$CHUNKS_TSV" | awk 'END { print NR + 0 }'
}

chunk_total_read_pairs() {
  if [[ ! -s "${CHUNKS_TSV:-}" ]]; then
    printf 'unknown\n'
    return
  fi
  awk -F '\t' '
    NR == 1 { next }
    $4 !~ /^[0-9]+$/ { available = 0; exit }
    { total += $4; available = 1 }
    END {
      if (available) {
        print total + 0
      } else {
        print "unknown"
      }
    }
  ' "$CHUNKS_TSV"
}

completed_chunk_count() {
  if [[ ! -s "${CHUNKS_TSV:-}" ]]; then
    printf '0\n'
    return
  fi
  local done_count=0
  local chunk_id chunk_r1 chunk_r2 n_read_pairs
  while IFS=$'\t' read -r chunk_id chunk_r1 chunk_r2 n_read_pairs; do
    [[ "$chunk_id" == "chunk_id" ]] && continue
    [[ -n "$chunk_id" ]] || continue
    if [[ "$(get_status "$(chunk_status_file "$chunk_id" chunk)")" == "done" ]]; then
      done_count=$((done_count + 1))
    fi
  done < "$CHUNKS_TSV"
  printf '%s\n' "$done_count"
}

print_run_summary() {
  progress_msg "run summary:"
  progress_msg "  sample: ${SAMPLE}"
  progress_msg "  input R1: ${R1}"
  progress_msg "  input R2: ${R2}"
  progress_msg "  genome name: ${GENOME_NAME}"
  progress_msg "  enzyme: ${ENZYME}"
  progress_msg "  threads: ${THREADS}"
  progress_msg "  chunk size: ${CHUNK_SIZE}"
  progress_msg "  max chunks: ${MAX_CHUNKS}"
  progress_msg "  workdir: ${WORKDIR}"
  progress_msg "  outdir: ${OUTDIR}"
}

print_chunk_split_summary() {
  local total_chunks total_pairs
  total_chunks="$(chunk_total_count)"
  total_pairs="$(chunk_total_read_pairs)"
  progress_msg "chunk splitting complete:"
  progress_msg "  total chunks: ${total_chunks}"
  progress_msg "  total read pairs: ${total_pairs}"
  progress_msg "  chunk size: ${CHUNK_SIZE}"
}

print_max_chunks_summary() {
  local processed_this_run="$1"
  local total_chunks completed_chunks
  total_chunks="$(chunk_total_count)"
  completed_chunks="$(completed_chunk_count)"
  progress_msg "processed ${processed_this_run} chunk(s) this run; ${completed_chunks}/${total_chunks} total chunks complete"
}

print_final_summary() {
  [[ -n "${STATUS_DIR:-}" ]] || return 0
  local total_chunks completed_chunks incomplete_chunks pipeline_status
  total_chunks="$(chunk_total_count)"
  completed_chunks="$(completed_chunk_count)"
  incomplete_chunks=$((total_chunks - completed_chunks))
  pipeline_status="$(get_status "${PIPELINE_STATUS:-}")"
  progress_msg "final summary:"
  progress_msg "  total chunks: ${total_chunks}"
  progress_msg "  completed chunks: ${completed_chunks}"
  progress_msg "  incomplete chunks: ${incomplete_chunks}"
  progress_msg "  final pipeline.status: ${pipeline_status}"
  progress_msg "  status directory: ${STATUS_DIR}"
  progress_msg "  logs directory: ${LOG_DIR}"
}

split_fastq_chunks() {
  if [[ "$(get_status "$CHUNK_SPLIT_STATUS")" == "done" ]] && validate_chunks_manifest; then
    log_msg "FASTQ chunks already complete"
    init_chunk_statuses
    print_chunk_split_summary
    return
  fi
  log_msg "Splitting paired FASTQ files into chunks"
  rm -rf "$CHUNK_FASTQ_DIR" "$CHUNK_STATUS_ROOT" "$CHUNKS_TSV"
  mkdir -p "$CHUNK_FASTQ_DIR" "$CHUNK_STATUS_ROOT"
  atomic_set_status "$CHUNK_SPLIT_STATUS" "running"
  CURRENT_STATUS_FILE="$CHUNK_SPLIT_STATUS"
  local log_file="${LOG_DIR}/split.log"
  local helper_manifest="${CHUNK_FASTQ_DIR}/chunks.tsv"
  local cmd_text="python3 ${SPLIT_HELPER} --r1 ${R1} --r2 ${R2} --chunk-size ${CHUNK_SIZE} --outdir ${CHUNK_FASTQ_DIR} --prefix chunk"
  write_step_log_header "$log_file" "$cmd_text" "$CHUNKS_TSV"
  set +e
  python3 "$SPLIT_HELPER" --r1 "$R1" --r2 "$R2" --chunk-size "$CHUNK_SIZE" --outdir "$CHUNK_FASTQ_DIR" --prefix chunk >> "$log_file" 2>&1
  local exit_code=$?
  set -e
  if [[ "$exit_code" -eq 0 && -s "$helper_manifest" ]]; then
    mv "$helper_manifest" "$CHUNKS_TSV"
  fi
  write_step_log_footer "$log_file" "$exit_code"
  CURRENT_STATUS_FILE=""
  if [[ "$exit_code" -eq 0 ]] && validate_chunks_manifest; then
    init_chunk_statuses
    atomic_set_status "$CHUNK_SPLIT_STATUS" "done"
    print_chunk_split_summary
  else
    atomic_set_status "$CHUNK_SPLIT_STATUS" "failed"
    die "FASTQ chunk splitting failed"
  fi
}

chunk_paths() {
  local chunk_id="$1"
  CHUNK_PROCESSED_DIR="${CHUNK_PROCESSED_ROOT}/${chunk_id}"
  CHUNK_WORK_DIR="${SAMPLE_WORK_CHUNKS}/${chunk_id}"
  TRIM_DIR="${CHUNK_PROCESSED_DIR}/trim"
  TRIM_R1="${TRIM_DIR}/${chunk_id}_val_1.fq.gz"
  TRIM_R2="${TRIM_DIR}/${chunk_id}_val_2.fq.gz"
  SAM_FILE="${CHUNK_WORK_DIR}/${chunk_id}.sam"
  SAM_TMP="${CHUNK_WORK_DIR}/${chunk_id}.sam.tmp"
  PARSED_FILE="${CHUNK_WORK_DIR}/${chunk_id}.parsed.pairsam.gz"
  SORTED_FILE="${CHUNK_WORK_DIR}/${chunk_id}.sorted.pairsam.gz"
  RESTRICTED_FILE="${CHUNK_WORK_DIR}/${chunk_id}.restricted.pairsam.gz"
  SELECTED_FILE="${CHUNK_PROCESSED_DIR}/chunk.selected.sorted.pairsam.gz"
  mkdir -p "$CHUNK_PROCESSED_DIR" "$CHUNK_WORK_DIR" "$TRIM_DIR"
}

chunk_status_file() {
  local chunk_id="$1"
  local step="$2"
  printf '%s/%s/%s.status\n' "$CHUNK_STATUS_ROOT" "$chunk_id" "$step"
}

delete_outputs_from_step() {
  local step="$1"
  case "$step" in
    trim_galore)
      rm -f "$TRIM_R1" "$TRIM_R2" "$SAM_FILE" "$SAM_TMP" "$PARSED_FILE" "$SORTED_FILE" "$RESTRICTED_FILE" "$SELECTED_FILE"
      ;;
    bwa_mem)
      rm -f "$SAM_FILE" "$SAM_TMP" "$PARSED_FILE" "$SORTED_FILE" "$RESTRICTED_FILE" "$SELECTED_FILE"
      ;;
    parse)
      rm -f "$PARSED_FILE" "$SORTED_FILE" "$RESTRICTED_FILE" "$SELECTED_FILE"
      ;;
    sort)
      rm -f "$SORTED_FILE" "$RESTRICTED_FILE" "$SELECTED_FILE"
      ;;
    restrict)
      rm -f "$RESTRICTED_FILE" "$SELECTED_FILE"
      ;;
    select)
      rm -f "$SELECTED_FILE"
      ;;
  esac
}

reset_step_and_downstream() {
  local chunk_id="$1"
  local step="$2"
  delete_outputs_from_step "$step"
  case "$step" in
    trim_galore)
      atomic_set_status "$(chunk_status_file "$chunk_id" trim_galore)" "null"
      atomic_set_status "$(chunk_status_file "$chunk_id" bwa_mem)" "null"
      atomic_set_status "$(chunk_status_file "$chunk_id" parse)" "null"
      atomic_set_status "$(chunk_status_file "$chunk_id" sort)" "null"
      atomic_set_status "$(chunk_status_file "$chunk_id" restrict)" "null"
      atomic_set_status "$(chunk_status_file "$chunk_id" select)" "null"
      ;;
    bwa_mem)
      atomic_set_status "$(chunk_status_file "$chunk_id" bwa_mem)" "null"
      atomic_set_status "$(chunk_status_file "$chunk_id" parse)" "null"
      atomic_set_status "$(chunk_status_file "$chunk_id" sort)" "null"
      atomic_set_status "$(chunk_status_file "$chunk_id" restrict)" "null"
      atomic_set_status "$(chunk_status_file "$chunk_id" select)" "null"
      ;;
    parse)
      atomic_set_status "$(chunk_status_file "$chunk_id" parse)" "null"
      atomic_set_status "$(chunk_status_file "$chunk_id" sort)" "null"
      atomic_set_status "$(chunk_status_file "$chunk_id" restrict)" "null"
      atomic_set_status "$(chunk_status_file "$chunk_id" select)" "null"
      ;;
    sort)
      atomic_set_status "$(chunk_status_file "$chunk_id" sort)" "null"
      atomic_set_status "$(chunk_status_file "$chunk_id" restrict)" "null"
      atomic_set_status "$(chunk_status_file "$chunk_id" select)" "null"
      ;;
    restrict)
      atomic_set_status "$(chunk_status_file "$chunk_id" restrict)" "null"
      atomic_set_status "$(chunk_status_file "$chunk_id" select)" "null"
      ;;
    select)
      atomic_set_status "$(chunk_status_file "$chunk_id" select)" "null"
      ;;
  esac
  atomic_set_status "$(chunk_status_file "$chunk_id" chunk)" "null"
}

validate_step_output() {
  local step="$1"
  case "$step" in
    trim_galore) validate_fastq_pair "$TRIM_R1" "$TRIM_R2" ;;
    bwa_mem) validate_sam "$SAM_FILE" ;;
    parse) validate_pairsam_gz "$PARSED_FILE" ;;
    sort) validate_pairsam_gz "$SORTED_FILE" ;;
    restrict) validate_pairsam_gz "$RESTRICTED_FILE" ;;
    select) validate_pairsam_gz "$SELECTED_FILE" ;;
    *) return 1 ;;
  esac
}

step_complete_or_reset() {
  local chunk_id="$1"
  local step="$2"
  local status_file
  status_file="$(chunk_status_file "$chunk_id" "$step")"
  if [[ "$(get_status "$status_file")" == "done" ]]; then
    if validate_step_output "$step"; then
      return 0
    fi
    log_msg "${chunk_id}: ${step} status was done but validation failed; resetting downstream"
    reset_step_and_downstream "$chunk_id" "$step"
    return 1
  fi
  reset_step_and_downstream "$chunk_id" "$step"
  return 1
}

run_trim_galore_step() {
  local chunk_id="$1"
  local chunk_r1="$2"
  local chunk_r2="$3"
  step_complete_or_reset "$chunk_id" "trim_galore" && return
  chunk_progress_msg "$chunk_id" "trim_galore"
  local status_file log_file cmd_text
  status_file="$(chunk_status_file "$chunk_id" trim_galore)"
  log_file="${CHUNK_LOG_DIR}/${chunk_id}.trim_galore.log"
  cmd_text="trim_galore -j ${THREADS} --paired --basename ${chunk_id} -o ${TRIM_DIR} ${chunk_r1} ${chunk_r2}"
  atomic_set_status "$status_file" "running"
  CURRENT_STATUS_FILE="$status_file"
  write_step_log_header "$log_file" "$cmd_text" "$TRIM_R1 $TRIM_R2"
  set +e
  trim_galore -j "$THREADS" --paired --basename "$chunk_id" -o "$TRIM_DIR" "$chunk_r1" "$chunk_r2" >> "$log_file" 2>&1
  local exit_code=$?
  set -e
  write_step_log_footer "$log_file" "$exit_code"
  CURRENT_STATUS_FILE=""
  if [[ "$exit_code" -eq 0 ]] && validate_fastq_pair "$TRIM_R1" "$TRIM_R2"; then
    atomic_set_status "$status_file" "done"
  else
    atomic_set_status "$status_file" "failed"
    die "${chunk_id}: trim_galore failed"
  fi
}

run_bwa_mem_step() {
  local chunk_id="$1"
  if [[ "$(get_status "$(chunk_status_file "$chunk_id" parse)")" == "done" && ! -s "$SAM_FILE" ]]; then
    return
  fi
  step_complete_or_reset "$chunk_id" "bwa_mem" && return
  chunk_progress_msg "$chunk_id" "bwa_mem"
  local status_file log_file cmd_text
  status_file="$(chunk_status_file "$chunk_id" bwa_mem)"
  log_file="${CHUNK_LOG_DIR}/${chunk_id}.bwa_mem.log"
  cmd_text="bwa mem -5SP -T0 -t ${THREADS} ${GENOME_FA} ${TRIM_R1} ${TRIM_R2} > ${SAM_TMP}; mv ${SAM_TMP} ${SAM_FILE}"
  atomic_set_status "$status_file" "running"
  CURRENT_STATUS_FILE="$status_file"
  write_step_log_header "$log_file" "$cmd_text" "$SAM_FILE"
  set +e
  bwa mem -5SP -T0 -t "$THREADS" "$GENOME_FA" "$TRIM_R1" "$TRIM_R2" > "$SAM_TMP" 2>> "$log_file"
  local exit_code=$?
  if [[ "$exit_code" -eq 0 ]]; then
    mv "$SAM_TMP" "$SAM_FILE"
    exit_code=$?
  fi
  set -e
  write_step_log_footer "$log_file" "$exit_code"
  CURRENT_STATUS_FILE=""
  if [[ "$exit_code" -eq 0 ]] && validate_sam "$SAM_FILE"; then
    atomic_set_status "$status_file" "done"
  else
    atomic_set_status "$status_file" "failed"
    die "${chunk_id}: bwa mem failed"
  fi
}

run_parse_step() {
  local chunk_id="$1"
  if [[ "$(get_status "$(chunk_status_file "$chunk_id" sort)")" == "done" && ! -s "$PARSED_FILE" ]]; then
    rm -f "$SAM_FILE"
    return
  fi
  step_complete_or_reset "$chunk_id" "parse" && { rm -f "$SAM_FILE"; return; }
  if [[ ! -s "$SAM_FILE" ]]; then
    atomic_set_status "$(chunk_status_file "$chunk_id" bwa_mem)" "null"
    run_bwa_mem_step "$chunk_id"
  fi
  chunk_progress_msg "$chunk_id" "parse"
  local status_file log_file cmd_text stats_file
  status_file="$(chunk_status_file "$chunk_id" parse)"
  log_file="${CHUNK_LOG_DIR}/${chunk_id}.parse.log"
  stats_file="${CHUNK_PROCESSED_DIR}/${chunk_id}.parse.stats"
  cmd_text="pairtools parse --min-mapq ${MIN_MAPQ} --walks-policy 5unique --max-inter-align-gap 30 --nproc-in ${THREADS} --nproc-out ${THREADS} --chroms-path ${GENOME_DIR}/${GENOME_NAME}.chrom.sizes --output-stats ${stats_file} --output ${PARSED_FILE} ${SAM_FILE}"
  atomic_set_status "$status_file" "running"
  CURRENT_STATUS_FILE="$status_file"
  write_step_log_header "$log_file" "$cmd_text" "$PARSED_FILE"
  set +e
  pairtools parse --min-mapq "$MIN_MAPQ" --walks-policy 5unique --max-inter-align-gap 30 \
    --nproc-in "$THREADS" --nproc-out "$THREADS" \
    --chroms-path "${GENOME_DIR}/${GENOME_NAME}.chrom.sizes" \
    --output-stats "$stats_file" --output "$PARSED_FILE" "$SAM_FILE" >> "$log_file" 2>&1
  local exit_code=$?
  set -e
  write_step_log_footer "$log_file" "$exit_code"
  CURRENT_STATUS_FILE=""
  if [[ "$exit_code" -eq 0 ]] && validate_pairsam_gz "$PARSED_FILE"; then
    atomic_set_status "$status_file" "done"
    rm -f "$SAM_FILE"
  else
    atomic_set_status "$status_file" "failed"
    die "${chunk_id}: pairtools parse failed"
  fi
}

run_sort_step() {
  local chunk_id="$1"
  if [[ "$(get_status "$(chunk_status_file "$chunk_id" restrict)")" == "done" && ! -s "$SORTED_FILE" ]]; then
    rm -f "$PARSED_FILE"
    return
  fi
  step_complete_or_reset "$chunk_id" "sort" && { rm -f "$PARSED_FILE"; return; }
  chunk_progress_msg "$chunk_id" "sort"
  local status_file log_file cmd_text
  status_file="$(chunk_status_file "$chunk_id" sort)"
  log_file="${CHUNK_LOG_DIR}/${chunk_id}.sort.log"
  cmd_text="pairtools sort --nproc ${THREADS} --tmpdir ${SAMPLE_TMPDIR} --output ${SORTED_FILE} ${PARSED_FILE}"
  atomic_set_status "$status_file" "running"
  CURRENT_STATUS_FILE="$status_file"
  write_step_log_header "$log_file" "$cmd_text" "$SORTED_FILE"
  set +e
  pairtools sort --nproc "$THREADS" --tmpdir "$SAMPLE_TMPDIR" --output "$SORTED_FILE" "$PARSED_FILE" >> "$log_file" 2>&1
  local exit_code=$?
  set -e
  write_step_log_footer "$log_file" "$exit_code"
  CURRENT_STATUS_FILE=""
  if [[ "$exit_code" -eq 0 ]] && validate_pairsam_gz "$SORTED_FILE"; then
    atomic_set_status "$status_file" "done"
    rm -f "$PARSED_FILE"
  else
    atomic_set_status "$status_file" "failed"
    die "${chunk_id}: pairtools sort failed"
  fi
}

run_restrict_step() {
  local chunk_id="$1"
  if [[ "$(get_status "$(chunk_status_file "$chunk_id" select)")" == "done" && ! -s "$RESTRICTED_FILE" ]]; then
    rm -f "$SORTED_FILE"
    return
  fi
  step_complete_or_reset "$chunk_id" "restrict" && { rm -f "$SORTED_FILE"; return; }
  chunk_progress_msg "$chunk_id" "restrict"
  local status_file log_file cmd_text frags_bed
  status_file="$(chunk_status_file "$chunk_id" restrict)"
  log_file="${CHUNK_LOG_DIR}/${chunk_id}.restrict.log"
  frags_bed="${GENOME_DIR}/${GENOME_NAME}.${ENZYME}.frags.bed"
  cmd_text="pairtools restrict -f ${frags_bed} --nproc-in ${THREADS} --nproc-out ${THREADS} --output ${RESTRICTED_FILE} ${SORTED_FILE}"
  atomic_set_status "$status_file" "running"
  CURRENT_STATUS_FILE="$status_file"
  write_step_log_header "$log_file" "$cmd_text" "$RESTRICTED_FILE"
  set +e
  pairtools restrict -f "$frags_bed" --nproc-in "$THREADS" --nproc-out "$THREADS" --output "$RESTRICTED_FILE" "$SORTED_FILE" >> "$log_file" 2>&1
  local exit_code=$?
  set -e
  write_step_log_footer "$log_file" "$exit_code"
  CURRENT_STATUS_FILE=""
  if [[ "$exit_code" -eq 0 ]] && validate_pairsam_gz "$RESTRICTED_FILE"; then
    atomic_set_status "$status_file" "done"
    rm -f "$SORTED_FILE"
  else
    atomic_set_status "$status_file" "failed"
    die "${chunk_id}: pairtools restrict failed"
  fi
}

run_select_step() {
  local chunk_id="$1"
  step_complete_or_reset "$chunk_id" "select" && { rm -f "$RESTRICTED_FILE"; return; }
  chunk_progress_msg "$chunk_id" "select"
  local status_file log_file cmd_text expr
  status_file="$(chunk_status_file "$chunk_id" select)"
  log_file="${CHUNK_LOG_DIR}/${chunk_id}.select.log"
  expr="(pair_type == 'UU') and not (chrom1 == chrom2 and abs(pos2 - pos1) < ${SHORT_CIS_CUTOFF})"
  cmd_text="pairtools select \"${expr}\" --output ${SELECTED_FILE} ${RESTRICTED_FILE}"
  atomic_set_status "$status_file" "running"
  CURRENT_STATUS_FILE="$status_file"
  write_step_log_header "$log_file" "$cmd_text" "$SELECTED_FILE"
  set +e
  pairtools select "$expr" --output "$SELECTED_FILE" "$RESTRICTED_FILE" >> "$log_file" 2>&1
  local exit_code=$?
  set -e
  write_step_log_footer "$log_file" "$exit_code"
  CURRENT_STATUS_FILE=""
  if [[ "$exit_code" -eq 0 ]] && validate_pairsam_gz "$SELECTED_FILE"; then
    atomic_set_status "$status_file" "done"
    rm -f "$RESTRICTED_FILE"
  else
    atomic_set_status "$status_file" "failed"
    die "${chunk_id}: pairtools select failed"
  fi
}

process_chunk() {
  local chunk_id="$1"
  local chunk_r1="$2"
  local chunk_r2="$3"
  local chunk_index="$4"
  local total_chunks="$5"
  CURRENT_CHUNK_INDEX="$chunk_index"
  CURRENT_CHUNK_TOTAL="$total_chunks"
  chunk_paths "$chunk_id"
  local chunk_status
  chunk_status="$(chunk_status_file "$chunk_id" chunk)"
  if [[ "$(get_status "$chunk_status")" == "done" ]] && validate_pairsam_gz "$SELECTED_FILE"; then
    chunk_progress_msg "$chunk_id" "already complete, skipping"
    return
  fi
  atomic_set_status "$chunk_status" "null"
  log_msg "${chunk_id}: processing"
  run_trim_galore_step "$chunk_id" "$chunk_r1" "$chunk_r2"
  run_bwa_mem_step "$chunk_id"
  run_parse_step "$chunk_id"
  run_sort_step "$chunk_id"
  run_restrict_step "$chunk_id"
  run_select_step "$chunk_id"
  if [[ "$(get_status "$(chunk_status_file "$chunk_id" select)")" == "done" ]] && validate_pairsam_gz "$SELECTED_FILE"; then
    atomic_set_status "$chunk_status" "done"
    log_msg "${chunk_id}: complete"
    chunk_progress_msg "$chunk_id" "done"
  else
    atomic_set_status "$chunk_status" "null"
    die "${chunk_id}: final chunk validation failed"
  fi
}

process_all_chunks() {
  local processed_this_run=0
  local total_chunks chunk_index=0
  total_chunks="$(chunk_total_count)"
  local chunk_id chunk_r1 chunk_r2 n_read_pairs
  while IFS=$'\t' read -r chunk_id chunk_r1 chunk_r2 n_read_pairs; do
    [[ "$chunk_id" == "chunk_id" ]] && continue
    [[ -n "$chunk_id" ]] || continue
    chunk_index=$((chunk_index + 1))
    CURRENT_CHUNK_INDEX="$chunk_index"
    CURRENT_CHUNK_TOTAL="$total_chunks"
    chunk_paths "$chunk_id"
    if [[ "$(get_status "$(chunk_status_file "$chunk_id" chunk)")" == "done" ]] && validate_pairsam_gz "$SELECTED_FILE"; then
      chunk_progress_msg "$chunk_id" "already complete, skipping"
      continue
    fi
    if [[ "$MAX_CHUNKS" -gt 0 && "$processed_this_run" -ge "$MAX_CHUNKS" ]]; then
      log_msg "Reached --max-chunks ${MAX_CHUNKS}; exiting cleanly"
      print_max_chunks_summary "$processed_this_run"
      return
    fi
    process_chunk "$chunk_id" "$chunk_r1" "$chunk_r2" "$chunk_index" "$total_chunks"
    processed_this_run=$((processed_this_run + 1))
  done < "$CHUNKS_TSV"
  if [[ "$MAX_CHUNKS" -gt 0 ]]; then
    print_max_chunks_summary "$processed_this_run"
  fi
}

check_all_chunks_done() {
  local total=0
  local chunk_id chunk_r1 chunk_r2 n_read_pairs
  while IFS=$'\t' read -r chunk_id chunk_r1 chunk_r2 n_read_pairs; do
    [[ "$chunk_id" == "chunk_id" ]] && continue
    [[ -n "$chunk_id" ]] || continue
    total=$((total + 1))
    chunk_paths "$chunk_id"
    if [[ "$(get_status "$(chunk_status_file "$chunk_id" chunk)")" != "done" ]] || ! validate_pairsam_gz "$SELECTED_FILE"; then
      atomic_set_status "$ALL_CHUNKS_STATUS" "null"
      atomic_set_status "$FINAL_STATUS" "null"
      return 1
    fi
  done < "$CHUNKS_TSV"
  [[ "$total" -gt 0 ]] || return 1
  atomic_set_status "$ALL_CHUNKS_STATUS" "done"
  atomic_set_status "$FINAL_STATUS" "done"
  atomic_set_status "$PIPELINE_STATUS" "done"
  {
    printf 'sample_id\t%s\n' "$SAMPLE"
    printf 'date\t%s\n' "$(date '+%Y-%m-%d %H:%M:%S')"
    printf 'number_of_chunks\t%s\n' "$total"
    printf 'chunk_size\t%s\n' "$CHUNK_SIZE"
    printf 'genome_name\t%s\n' "$GENOME_NAME"
    printf 'enzyme\t%s\n' "$ENZYME"
    printf 'min_mapq\t%s\n' "$MIN_MAPQ"
    printf 'short_cis_cutoff\t%s\n' "$SHORT_CIS_CUTOFF"
  } > "${STATUS_DIR}/pipeline.done"
  return 0
}

print_status() {
  local total=0 done_count=0 incomplete_count=0
  local trim_incomplete=0 bwa_incomplete=0 parse_incomplete=0 sort_incomplete=0 restrict_incomplete=0 select_incomplete=0
  if [[ -s "$CHUNKS_TSV" ]]; then
    local chunk_id chunk_r1 chunk_r2 n_read_pairs
    while IFS=$'\t' read -r chunk_id chunk_r1 chunk_r2 n_read_pairs; do
      [[ "$chunk_id" == "chunk_id" ]] && continue
      [[ -n "$chunk_id" ]] || continue
      total=$((total + 1))
      if [[ "$(get_status "$(chunk_status_file "$chunk_id" chunk)")" == "done" ]]; then
        done_count=$((done_count + 1))
      else
        incomplete_count=$((incomplete_count + 1))
      fi
      [[ "$(get_status "$(chunk_status_file "$chunk_id" trim_galore)")" == "done" ]] || trim_incomplete=$((trim_incomplete + 1))
      [[ "$(get_status "$(chunk_status_file "$chunk_id" bwa_mem)")" == "done" ]] || bwa_incomplete=$((bwa_incomplete + 1))
      [[ "$(get_status "$(chunk_status_file "$chunk_id" parse)")" == "done" ]] || parse_incomplete=$((parse_incomplete + 1))
      [[ "$(get_status "$(chunk_status_file "$chunk_id" sort)")" == "done" ]] || sort_incomplete=$((sort_incomplete + 1))
      [[ "$(get_status "$(chunk_status_file "$chunk_id" restrict)")" == "done" ]] || restrict_incomplete=$((restrict_incomplete + 1))
      [[ "$(get_status "$(chunk_status_file "$chunk_id" select)")" == "done" ]] || select_incomplete=$((select_incomplete + 1))
    done < "$CHUNKS_TSV"
  fi
  local percent_complete="0.0"
  if [[ "$total" -gt 0 ]]; then
    percent_complete="$(awk -v done="$done_count" -v total="$total" 'BEGIN { printf "%.1f", (done / total) * 100 }')"
  fi
  cat <<EOF
Sample: ${SAMPLE}
Pipeline: $(get_status "$PIPELINE_STATUS")
QC: $(get_status "$QC_STATUS")
Chunk split: $(get_status "$CHUNK_SPLIT_STATUS")
All chunks: $(get_status "$ALL_CHUNKS_STATUS")
Final: $(get_status "$FINAL_STATUS")

Chunks:
  total: ${total}
  done: ${done_count}
  incomplete: ${incomplete_count}
  complete: ${done_count}/${total} (${percent_complete}%)

Per-step incomplete counts:
  trim_galore: ${trim_incomplete}
  bwa_mem: ${bwa_incomplete}
  parse: ${parse_incomplete}
  sort: ${sort_incomplete}
  restrict: ${restrict_incomplete}
  select: ${select_incomplete}
EOF
}

acquire_lock() {
  mkdir -p "$LOCK_ROOT"
  local lock_file="${LOCK_ROOT}/${SAMPLE}.lock"
  if command -v flock >/dev/null 2>&1; then
    eval "exec ${LOCK_FD}>\"${lock_file}\""
    flock -n "$LOCK_FD" || die "Another pipeline instance is already processing sample ${SAMPLE}"
    LOCK_HELD=1
  else
    LOCK_DIR="${lock_file}.d"
    if ! mkdir "$LOCK_DIR" 2>/dev/null; then
      die "Another pipeline instance may be processing sample ${SAMPLE}; lock exists: ${LOCK_DIR}"
    fi
    LOCK_HELD=1
  fi
}

acquire_reference_lock() {
  mkdir -p "${GENOME_DIR}/.locks"
  local lock_file="${GENOME_DIR}/.locks/${GENOME_NAME}.${ENZYME}.reference.lock"
  if command -v flock >/dev/null 2>&1; then
    eval "exec ${REF_LOCK_FD}>\"${lock_file}\""
    flock "$REF_LOCK_FD" || die "Could not acquire reference lock for ${GENOME_NAME}/${ENZYME}"
    REF_LOCK_HELD=1
  else
    REF_LOCK_DIR="${lock_file}.d"
    while ! mkdir "$REF_LOCK_DIR" 2>/dev/null; do
      log_msg "Waiting for reference lock: ${REF_LOCK_DIR}"
      sleep 10
    done
    REF_LOCK_HELD=1
  fi
}

release_reference_lock() {
  if [[ "$REF_LOCK_HELD" -ne 1 ]]; then
    return
  fi
  if command -v flock >/dev/null 2>&1; then
    flock -u "$REF_LOCK_FD" || true
  fi
  if [[ -n "${REF_LOCK_DIR:-}" && -d "$REF_LOCK_DIR" ]]; then
    rm -rf "$REF_LOCK_DIR"
    REF_LOCK_DIR=""
  fi
  REF_LOCK_HELD=0
}

cleanup_lock() {
  if [[ -n "${LOCK_DIR:-}" && -d "$LOCK_DIR" ]]; then
    rm -rf "$LOCK_DIR"
  fi
  LOCK_HELD=0
  release_reference_lock
}

preflight() {
  [[ -r "$R1" ]] || die "R1 FASTQ is missing or unreadable: ${R1}"
  [[ -r "$R2" ]] || die "R2 FASTQ is missing or unreadable: ${R2}"
  [[ -r "$GENOME_FA" ]] || die "Genome FASTA is missing or unreadable: ${GENOME_FA}"
  require_command fastqc
  require_command trim_galore
  require_command bwa
  require_command samtools
  require_command pairtools
  require_command cooler
  require_command gzip
  require_command python3
  [[ -s "$SPLIT_HELPER" ]] || die "FASTQ split helper not found: ${SPLIT_HELPER}"
}

parse_args() {
  while [[ $# -gt 0 ]]; do
    case "$1" in
      --sample) SAMPLE="$2"; shift 2 ;;
      --r1) R1="$2"; shift 2 ;;
      --r2) R2="$2"; shift 2 ;;
      --genome-name) GENOME_NAME="$2"; shift 2 ;;
      --genome-fa) GENOME_FA="$2"; shift 2 ;;
      --enzyme) ENZYME="$2"; shift 2 ;;
      --workdir) WORKDIR="$2"; shift 2 ;;
      --outdir) OUTDIR="$2"; shift 2 ;;
      --threads) THREADS="$2"; shift 2 ;;
      --chunk-size) CHUNK_SIZE="$2"; shift 2 ;;
      --min-mapq) MIN_MAPQ="$2"; shift 2 ;;
      --short-cis-cutoff) SHORT_CIS_CUTOFF="$2"; shift 2 ;;
      --max-chunks) MAX_CHUNKS="$2"; shift 2 ;;
      --force-init) FORCE_INIT=1; shift ;;
      --status-only) STATUS_ONLY=1; shift ;;
      -h|--help) usage; exit 0 ;;
      *) usage >&2; exit 2 ;;
    esac
  done
  [[ -n "$SAMPLE" && -n "$R1" && -n "$R2" && -n "$GENOME_NAME" && -n "$GENOME_FA" && -n "$ENZYME" && -n "$WORKDIR" && -n "$OUTDIR" ]] || {
    usage >&2
    exit 2
  }
  [[ "$THREADS" =~ ^[0-9]+$ && "$THREADS" -gt 0 ]] || die "--threads must be a positive integer"
  [[ "$CHUNK_SIZE" =~ ^[0-9]+$ && "$CHUNK_SIZE" -gt 0 ]] || die "--chunk-size must be a positive integer"
  [[ "$MIN_MAPQ" =~ ^[0-9]+$ ]] || die "--min-mapq must be a nonnegative integer"
  [[ "$SHORT_CIS_CUTOFF" =~ ^[0-9]+$ ]] || die "--short-cis-cutoff must be a nonnegative integer"
  [[ "$MAX_CHUNKS" =~ ^[0-9]+$ ]] || die "--max-chunks must be a nonnegative integer"
}

main() {
  parse_args "$@"
  prepare_dirs
  trap on_error ERR
  trap on_interrupt INT TERM
  trap cleanup_lock EXIT
  if [[ "$STATUS_ONLY" -eq 1 ]]; then
    print_status
    exit 0
  fi
  acquire_lock
  init_status
  print_run_summary
  preflight
  atomic_set_status "$PIPELINE_STATUS" "running"
  prepare_reference
  run_fastqc
  split_fastq_chunks
  process_all_chunks
  if check_all_chunks_done; then
    log_msg "Pipeline complete for sample ${SAMPLE}"
  else
    atomic_set_status "$PIPELINE_STATUS" "null"
    log_msg "Pipeline stopped cleanly with incomplete chunks; pipeline status reset to null"
  fi
  print_final_summary
}

main "$@"
