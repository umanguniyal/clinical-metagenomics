#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG_YAML="${SCRIPT_DIR}/config/config.yaml"

CORES=8
INTERACTIVE_MANIFEST="false"
ASK_PROFILE="true"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --cores) CORES="${2:-8}"; shift 2 ;;
    --interactive-manifest) INTERACTIVE_MANIFEST="true"; shift ;;
    --no-ask-profile) ASK_PROFILE="false"; shift ;;
    *) echo "Unknown arg: $1"; echo "Usage: ./pipeline run [--cores N] [--interactive-manifest] [--no-ask-profile]"; exit 1 ;;
  esac
done

PROJECT_ROOT="$(python3 -c 'import os,yaml; cfg=yaml.safe_load(open("'"$CONFIG_YAML"'")); root=cfg.get("project_root",".."); print(os.path.realpath(os.path.join(os.getcwd(), root)))')"
MANIFEST="${PROJECT_ROOT}/manifests/samples.tsv"

mkdir -p "${PROJECT_ROOT}/manifests"

if [[ "${INTERACTIVE_MANIFEST}" == "true" ]]; then
  python3 "${SCRIPT_DIR}/scripts/manifest_cli.py" --out "${MANIFEST}"
fi

if [[ ! -f "${MANIFEST}" ]]; then
  echo "ERROR: Manifest not found: ${MANIFEST}"
  echo "Run: ./pipeline ingest"
  exit 1
fi

echo "Using manifest: ${MANIFEST}"
echo "Manifest preview:"
column -t -s $'\t' "${MANIFEST}" | sed -n '1,40p'
echo ""

PROFILE_OVERRIDE=""
if [[ "${ASK_PROFILE}" == "true" ]]; then
  echo "=== PROFILE ==="
  echo "  1) bioinformatician"
  echo "  2) clinician"
  read -rp "Choice [1/2] (default: 1): " PR
  PR="${PR:-1}"
  case "${PR}" in
    1) PROFILE_OVERRIDE="bioinformatician" ;;
    2) PROFILE_OVERRIDE="clinician" ;;
    *) echo "Invalid profile choice"; exit 1 ;;
  esac
  echo "Selected profile: ${PROFILE_OVERRIDE}"
  echo ""
fi

cd "${SCRIPT_DIR}"

if [[ -n "${PROFILE_OVERRIDE}" ]]; then
  snakemake --use-conda --rerun-incomplete --cores "${CORES}" all --config profile="${PROFILE_OVERRIDE}"
else
  snakemake --use-conda --rerun-incomplete --cores "${CORES}" all
fi
