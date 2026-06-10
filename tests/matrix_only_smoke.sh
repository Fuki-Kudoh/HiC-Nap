#!/usr/bin/env bash
set -Eeuo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TMPDIR="$(mktemp -d)"
trap 'rm -rf "$TMPDIR"' EXIT

BIN_DIR="${TMPDIR}/bin"
OUTDIR="${TMPDIR}/out"
NO_HIC_OUTDIR="${TMPDIR}/out-no-hic"
GENOME_FA="${TMPDIR}/tiny.fa"
mkdir -p "$BIN_DIR" "${OUTDIR}/pairs" "${NO_HIC_OUTDIR}/pairs"

cat > "${BIN_DIR}/pairix" <<'STUB'
#!/usr/bin/env bash
set -Eeuo pipefail
if [[ "${1:-}" == "-l" ]]; then
  printf 'chrTiny\n'
  exit 0
fi
pairs_path="${1:?pairix needs a pairs file}"
: > "${pairs_path}.px2"
STUB

cat > "${BIN_DIR}/cooler" <<'STUB'
#!/usr/bin/env bash
set -Eeuo pipefail
case "${1:-}" in
  info|ls)
    [[ -s "${2:?cooler validation needs a file}" ]]
    ;;
  cload)
    out="${@: -1}"
    printf 'cool\n' > "$out"
    ;;
  zoomify)
    out=""
    while [[ $# -gt 0 ]]; do
      if [[ "$1" == "-o" ]]; then
        out="$2"
        break
      fi
      shift
    done
    [[ -n "$out" ]]
    printf 'mcool\n' > "$out"
    ;;
  *)
    printf 'unexpected cooler command: %s\n' "$*" >&2
    exit 1
    ;;
esac
STUB

cat > "${BIN_DIR}/juicer_tools" <<'STUB'
#!/usr/bin/env bash
set -Eeuo pipefail
[[ "${1:-}" == "pre" ]] || exit 1
shift
if [[ "${1:-}" == "-j" ]]; then
  shift 2
fi
valid_pairs="${1:?missing valid pairs}"
hic_file="${2:?missing hic output}"
chrom_sizes="${3:?missing chrom sizes}"
[[ -s "$valid_pairs" && -s "$chrom_sizes" ]]
printf 'hic\n' > "$hic_file"
STUB

cat > "${BIN_DIR}/samtools" <<'STUB'
#!/usr/bin/env bash
set -Eeuo pipefail
if [[ "${1:-}" == "faidx" ]]; then
  printf 'chrTiny\t1000\t9\t80\t81\n' > "${2}.fai"
fi
STUB

chmod +x "${BIN_DIR}/pairix" "${BIN_DIR}/cooler" "${BIN_DIR}/juicer_tools" "${BIN_DIR}/samtools"

cat > "$GENOME_FA" <<'EOF_FASTA'
>chrTiny
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
EOF_FASTA
printf 'chrTiny\t80\t9\t80\t81\n' > "${GENOME_FA}.fai"

cat > "${TMPDIR}/tiny.valid.pairs" <<'EOF_PAIRS'
## pairs format v1.0
#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type
read1 chrTiny 10 chrTiny 40 + - UU
EOF_PAIRS
gzip -c "${TMPDIR}/tiny.valid.pairs" > "${OUTDIR}/pairs/TINY.valid.pairs.gz"
printf 'pairix-index\n' > "${OUTDIR}/pairs/TINY.valid.pairs.gz.px2"
cp "${OUTDIR}/pairs/TINY.valid.pairs.gz" "${NO_HIC_OUTDIR}/pairs/TINY.valid.pairs.gz"
cp "${OUTDIR}/pairs/TINY.valid.pairs.gz.px2" "${NO_HIC_OUTDIR}/pairs/TINY.valid.pairs.gz.px2"

PATH="${BIN_DIR}:${PATH}" "${REPO_ROOT}/hic_chunk_pipeline.sh" \
  --sample TINY \
  --genome-name tiny \
  --genome-fa "$GENOME_FA" \
  --outdir "$OUTDIR" \
  --threads 2 \
  --matrix-only

[[ -s "${OUTDIR}/cool/TINY.10000.cool" ]]
[[ -s "${OUTDIR}/cool/TINY.mcool" ]]
[[ -s "${OUTDIR}/hic/TINY.hic" ]]
grep -q '^done$' "${OUTDIR}/status/TINY/matrix.status"
grep -q '^cool	' "${OUTDIR}/status/TINY/pipeline.done"
grep -q '^mcool	' "${OUTDIR}/status/TINY/pipeline.done"
grep -q '^hic	' "${OUTDIR}/status/TINY/pipeline.done"

PATH="${BIN_DIR}:${PATH}" "${REPO_ROOT}/hic_chunk_pipeline.sh" \
  --sample TINY \
  --genome-name tiny \
  --genome-fa "$GENOME_FA" \
  --outdir "$NO_HIC_OUTDIR" \
  --threads 2 \
  --matrix-only \
  --no-hic

[[ -s "${NO_HIC_OUTDIR}/cool/TINY.10000.cool" ]]
[[ -s "${NO_HIC_OUTDIR}/cool/TINY.mcool" ]]
[[ ! -e "${NO_HIC_OUTDIR}/hic/TINY.hic" ]]
grep -q '^done$' "${NO_HIC_OUTDIR}/status/TINY/matrix.status"
grep -q '^cool	' "${NO_HIC_OUTDIR}/status/TINY/pipeline.done"
grep -q '^mcool	' "${NO_HIC_OUTDIR}/status/TINY/pipeline.done"
! grep -q '^hic	' "${NO_HIC_OUTDIR}/status/TINY/pipeline.done"
