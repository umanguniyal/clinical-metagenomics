#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG_YAML="${SCRIPT_DIR}/config/config.yaml"

CORES=8
AUTO_RUN="false"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --run) AUTO_RUN="true"; shift ;;
    --cores) CORES="${2:-8}"; shift 2 ;;
    *) echo "Unknown arg: $1"; echo "Usage: ./run_ingest.sh [--run] [--cores N]"; exit 1 ;;
  esac
done

PROJECT_ROOT="$(python3 -c 'import os,yaml; cfg=yaml.safe_load(open("'"$CONFIG_YAML"'")); print(os.path.realpath(os.path.join(os.getcwd(), cfg.get("project_root",".."))))')"

MANIFEST_DIR="${PROJECT_ROOT}/manifests"
MANIFEST="${MANIFEST_DIR}/samples.tsv"

mkdir -p "${MANIFEST_DIR}" "${PROJECT_ROOT}/illumina" "${PROJECT_ROOT}/nanopore"

echo ""
echo "Project root: ${PROJECT_ROOT}"
echo "Manifest will be written to: ${MANIFEST}"
echo ""

echo "=== PLATFORM (single-platform per run) ==="
echo "  1) illumina (paired-end)"
echo "  2) nanopore (single-end)"
read -rp "Choice [1/2] (default: 1): " PCH
PCH="${PCH:-1}"

case "${PCH}" in
  1) RUN_PLATFORM="illumina" ;;
  2) RUN_PLATFORM="nanopore" ;;
  *) echo "Invalid platform choice"; exit 1 ;;
esac

echo "Setting config platform to: ${RUN_PLATFORM}"
python3 - <<PY
import yaml
p = "${CONFIG_YAML}"
cfg = yaml.safe_load(open(p))
cfg["platform"] = "${RUN_PLATFORM}"
with open(p, "w") as f:
    yaml.safe_dump(cfg, f, sort_keys=False)
PY

echo ""
echo "=== INGEST OPTIONS ==="
echo "1) Paste URL(s) (one per line; blank line to finish)"
echo "2) Provide a link file path (one URL/accession/path per line)"
echo "3) Paste accession(s) (ERR/SRR/DRR/SRX/ERX/DRX; one per line; blank line to finish)"
echo "4) Local path (file or directory) containing FASTQ/FASTQ.GZ"
echo "5) Skip ingest: only scan existing illumina/ and nanopore/ and write manifest"
echo ""
read -rp "Choice [1/2/3/4/5]: " CH

INPUT_LIST=()

if [[ "${CH}" == "1" ]]; then
  echo "Paste URLs (blank line ends):"
  while IFS= read -r line; do
    [[ -z "${line}" ]] && break
    INPUT_LIST+=("${line}")
  done
elif [[ "${CH}" == "2" ]]; then
  read -rp "Path to link file: " LF
  LF="${LF/#\~/$HOME}"
  [[ ! -f "${LF}" ]] && { echo "Link file not found: ${LF}"; exit 1; }
  mapfile -t INPUT_LIST < <(grep -vE '^\s*$|^\s*#' "${LF}")
elif [[ "${CH}" == "3" ]]; then
  echo "Paste accessions (blank line ends):"
  while IFS= read -r line; do
    [[ -z "${line}" ]] && break
    INPUT_LIST+=("${line}")
  done
elif [[ "${CH}" == "4" ]]; then
  read -rp "Local file/dir path: " P
  P="${P/#\~/$HOME}"
  INPUT_LIST+=("${P}")
elif [[ "${CH}" == "5" ]]; then
  INPUT_LIST=()
else
  echo "Invalid choice."
  exit 1
fi

# Only pass --inputs if there are actual inputs
if [[ ${#INPUT_LIST[@]} -gt 0 ]]; then
  python3 "${SCRIPT_DIR}/scripts/ingest_stage.py" \
    --project-root "${PROJECT_ROOT}" \
    --out-manifest "${MANIFEST}" \
    --platform "${RUN_PLATFORM}" \
    --inputs "${INPUT_LIST[@]}"
else
  python3 "${SCRIPT_DIR}/scripts/ingest_stage.py" \
    --project-root "${PROJECT_ROOT}" \
    --out-manifest "${MANIFEST}" \
    --platform "${RUN_PLATFORM}"
fi

echo ""
echo "DONE. Manifest written: ${MANIFEST}"
echo "Preview (first 200 lines):"
column -t -s $'\t' "${MANIFEST}" | sed -n '1,200p'

sync || true
echo ""
echo "Verifying manifest exists and has at least 2 lines (header + >=1 sample)..."

for i in {1..80}; do
  if [[ -f "${MANIFEST}" ]]; then
    LINES="$(wc -l < "${MANIFEST}" | tr -d ' ')"
    if [[ "${LINES}" -ge 2 ]]; then
      echo "Manifest OK (lines=${LINES})."
      break
    fi
  fi
  sleep 0.25
done

if [[ ! -f "${MANIFEST}" ]]; then
  echo "ERROR: Manifest not found after ingest: ${MANIFEST}"
  exit 1
fi

LINES="$(wc -l < "${MANIFEST}" | tr -d ' ')"
if [[ "${LINES}" -lt 2 ]]; then
  echo "ERROR: Manifest exists but appears empty (lines=${LINES}): ${MANIFEST}"
  exit 1
fi

if [[ "${AUTO_RUN}" != "true" ]]; then
  echo ""
  read -rp "Run Snakemake now? [y/N]: " GO
  [[ "${GO}" =~ ^[Yy]$ ]] && AUTO_RUN="true"
fi

if [[ "${AUTO_RUN}" == "true" ]]; then
  echo ""
  echo "Running pipeline via run_smk.sh with cores=${CORES}"
  exec "${SCRIPT_DIR}/run_smk.sh" --cores "${CORES}"
fi
