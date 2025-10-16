#!/bin/bash
# Upload methyltna references to Google Cloud Storage
# Run this script after authenticating with: gcloud auth login

set -e

echo "=== Uploading methyltna references to GCS ==="
echo ""

GCS_BASE="gs://motleybio/Workspaces/Marcus/methyltna_references"

# Check if we can access the bucket
echo "Testing GCS access..."
if ! gsutil ls "${GCS_BASE}/" &>/dev/null; then
    echo "Creating base directory..."
    echo "" | gsutil cp - "${GCS_BASE}/.placeholder" 2>/dev/null || true
fi

echo ""
echo "Starting uploads (this will take a while with 55GB of data)..."
echo ""

# Upload README (4KB - quick)
echo "[1/6] Uploading README.md..."
gsutil -m cp /home/marcus/pipelines/methyltna/references/README.md "${GCS_BASE}/"

# Upload GTF (1.5GB - ~1 minute)
echo "[2/6] Uploading GTF annotation (1.5GB)..."
gsutil -m cp /home/marcus/pipelines/methyltna/references/Homo_sapiens.GRCh38.112.chr_label.gtf "${GCS_BASE}/"

# Upload FASTA (3.5GB - ~2 minutes)
echo "[3/6] Uploading FASTA files (3.5GB)..."
gsutil -m cp -r /home/marcus/pipelines/methyltna/references/fasta "${GCS_BASE}/"

# Upload Bowtie2 indexes (8.7GB - ~5 minutes)
echo "[4/6] Uploading Bowtie2 indexes (8.7GB)..."
gsutil -m cp -r /home/marcus/pipelines/methyltna/references/bowtie2_indexes "${GCS_BASE}/"

# Upload Biscuit indexes (9.8GB - ~5 minutes)
echo "[5/6] Uploading Biscuit indexes (9.8GB)..."
gsutil -m cp -r /home/marcus/pipelines/methyltna/references/biscuit_indexes "${GCS_BASE}/"

# Upload STAR indexes (30GB - ~15 minutes)
echo "[6/6] Uploading STAR indexes (30GB)..."
gsutil -m cp -r /home/marcus/pipelines/methyltna/references/star_indexes "${GCS_BASE}/"

echo ""
echo "=== Upload Complete! ==="
echo ""
echo "Uploaded references to: ${GCS_BASE}"
echo ""
echo "Verifying uploads..."
gsutil ls -lh "${GCS_BASE}/" | grep -E "TOTAL|star_indexes|biscuit_indexes|bowtie2_indexes|fasta|gtf|README"
echo ""
echo "✅ All references uploaded successfully!"
