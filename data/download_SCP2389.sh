#!/usr/bin/env bash
# Download script for Broad Institute Single Cell Portal SCP2389
set -euo pipefail

TARGET_DIR="data/SCP2389/expression"
mkdir -p "${TARGET_DIR}"

echo "Downloading Seq-Well raw expression files for SCP2389..."
# Sostituire con il link curl generato dalla sessione del Single Cell Portal:
# curl -k -o "${TARGET_DIR}/matrix.mtx" "<URL_DOWNLOAD_MTX>"
# curl -k -o "${TARGET_DIR}/barcodes.tsv" "<URL_DOWNLOAD_BARCODES>"
# curl -k -o "${TARGET_DIR}/genes.tsv" "<URL_DOWNLOAD_GENES>"

echo "Download completed."
