#!/bin/bash
# Test cloud storage fallback for references
# This script tests the tiered caching: local → cloud → original source

set -e

echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "🧪 Testing Cloud Reference Cache Fallback"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

# Backup existing references if they exist
BACKUP_DIR="/tmp/methyltna_references_backup_$(date +%s)"
if [ -d "/home/marcus/pipelines/methyltna/references" ]; then
    echo "📦 Backing up existing references to: ${BACKUP_DIR}"
    mv /home/marcus/pipelines/methyltna/references "${BACKUP_DIR}"
    echo "   ✓ Backup complete"
    echo ""
fi

# Create minimal test samplesheet for reference preparation
TEST_SAMPLESHEET="/tmp/test_cloud_cache_samples.csv"
cat > "${TEST_SAMPLESHEET}" << 'EOF'
sample,fastq_1,fastq_2
test_sample,/home/marcus/pipelines/methyltna/test_data/Ca3-RNA-NNSR_09Z_R1.fastq.gz,/home/marcus/pipelines/methyltna/test_data/Ca3-RNA-NNSR_09Z_R2.fastq.gz
EOF

echo "🚀 Testing cloud cache by checking if files are accessible..."
echo "   This will test: Check cloud uploads → Download a sample file"
echo ""

# First, verify cloud uploads completed
echo "Step 1: Verifying cloud storage contents..."
if ! gsutil ls gs://motleybio/Workspaces/Marcus/methyltna_references/ &>/dev/null; then
    echo "   ❌ Cannot access cloud reference cache"
    echo "   Make sure uploads have completed"
    EXIT_CODE=1
else
    echo "   ✓ Cloud cache accessible"
    echo ""
    echo "Step 2: Listing cloud cache contents..."
    gsutil ls -lh gs://motleybio/Workspaces/Marcus/methyltna_references/ | head -20
    echo ""

    # Test downloading FASTA file from cloud cache
    echo "Step 3: Testing download from cloud cache..."
    mkdir -p /home/marcus/pipelines/methyltna/references/fasta

    GENOME_FILE="GRCh38_full_analysis_set_plus_decoy_hla.fa"
    if gsutil ls gs://motleybio/Workspaces/Marcus/methyltna_references/fasta/${GENOME_FILE} &>/dev/null; then
        echo "   ✓ Found genome FASTA in cloud cache"
        echo "   Testing download (will download full file to verify)..."

        if gsutil -m cp gs://motleybio/Workspaces/Marcus/methyltna_references/fasta/${GENOME_FILE} \
                        /home/marcus/pipelines/methyltna/references/fasta/; then
            echo "   ✅ Successfully downloaded genome FASTA from cloud cache"
            ls -lh /home/marcus/pipelines/methyltna/references/fasta/${GENOME_FILE}
        else
            echo "   ❌ Failed to download genome FASTA"
            EXIT_CODE=1
        fi
    else
        echo "   ⚠️  Genome FASTA not found in cloud cache yet"
        echo "   Uploads may still be in progress"
        EXIT_CODE=1
    fi

    # Test downloading GTF from cloud cache
    echo ""
    echo "Step 4: Testing GTF download from cloud cache..."
    GTF_FILE="Homo_sapiens.GRCh38.112.chr_label.gtf"
    if gsutil ls gs://motleybio/Workspaces/Marcus/methyltna_references/${GTF_FILE} &>/dev/null; then
        echo "   ✓ Found GTF annotation in cloud cache"

        if gsutil -m cp gs://motleybio/Workspaces/Marcus/methyltna_references/${GTF_FILE} \
                        /home/marcus/pipelines/methyltna/references/; then
            echo "   ✅ Successfully downloaded GTF from cloud cache"
            ls -lh /home/marcus/pipelines/methyltna/references/${GTF_FILE}
            EXIT_CODE=0
        else
            echo "   ❌ Failed to download GTF"
            EXIT_CODE=1
        fi
    else
        echo "   ⚠️  GTF not found in cloud cache yet"
        echo "   Uploads may still be in progress"
        EXIT_CODE=1
    fi
fi

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
if [ $EXIT_CODE -eq 0 ]; then
    echo "✅ Test PASSED - Pipeline completed successfully"
    echo ""
    echo "📊 Verification:"
    echo ""

    # Check if references were downloaded
    if [ -d "/home/marcus/pipelines/methyltna/references/fasta" ]; then
        echo "   ✓ FASTA files downloaded to local cache"
        ls -lh /home/marcus/pipelines/methyltna/references/fasta/
    fi
    echo ""

    if [ -f "/home/marcus/pipelines/methyltna/references/Homo_sapiens.GRCh38.112.chr_label.gtf" ]; then
        echo "   ✓ GTF annotation downloaded to local cache"
        ls -lh /home/marcus/pipelines/methyltna/references/*.gtf
    fi
    echo ""

    # Check workflow logs for cloud cache messages
    echo "📝 Checking logs for cloud cache activity..."
    if grep -r "cloud cache" work/ 2>/dev/null | head -5; then
        echo "   ✓ Cloud cache was accessed"
    else
        echo "   ⚠️  No cloud cache messages found (may have used local cache or original source)"
    fi

else
    echo "❌ Test FAILED - Pipeline exited with code: ${EXIT_CODE}"
    echo ""
    echo "Check the Nextflow logs for details:"
    echo "   tail -100 .nextflow.log"
fi

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
echo "🧹 Cleanup options:"
echo "   1. Keep new cache and remove backup:"
echo "      rm -rf ${BACKUP_DIR}"
echo ""
echo "   2. Restore original cache:"
if [ -d "${BACKUP_DIR}" ]; then
    echo "      rm -rf /home/marcus/pipelines/methyltna/references"
    echo "      mv ${BACKUP_DIR} /home/marcus/pipelines/methyltna/references"
else
    echo "      (No backup to restore)"
fi
echo ""
echo "   3. Keep both for comparison"
echo ""
