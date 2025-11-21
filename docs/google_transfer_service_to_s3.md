# Google Transfer Service: GCS to AWS S3 Migration Guide

**Purpose:** Transfer multi-TB pipeline data from Google Cloud Storage to AWS S3
**Use Case:** Migrating dozens of terabytes of pipeline outputs for methyltna project
**Cost Advantage:** **Zero GCS egress charges** (saves ~$120/TB compared to other methods)

---

## Table of Contents
- [Overview](#overview)
- [Cost Analysis](#cost-analysis)
- [Prerequisites](#prerequisites)
- [AWS Setup](#aws-setup)
- [GCP Setup](#gcp-setup)
- [Creating the Transfer Job](#creating-the-transfer-job)
- [Monitoring the Transfer](#monitoring-the-transfer)
- [Validation](#validation)
- [Troubleshooting](#troubleshooting)
- [Best Practices](#best-practices)

---

## Overview

**Google Cloud Storage Transfer Service to AWS S3** is Google's managed service for transferring data between GCS and S3. It's the most cost-effective method for large-scale migrations.

### Why Choose This Method?

For transferring **dozens of terabytes**, this is the clear winner:

| Transfer Size | Traditional Method Cost | Google Transfer Service Cost | Savings |
|---------------|------------------------|------------------------------|---------|
| 10 TB         | $1,200 (egress)        | $0                          | $1,200  |
| 50 TB         | $6,000 (egress)        | $0                          | $6,000  |
| 100 TB        | $12,000 (egress)       | $0                          | $12,000 |

**Additional Benefits:**
- ✅ **Zero egress charges** from GCS
- ✅ **Managed service** - Google handles retries, validation, progress tracking
- ✅ **Optimized network routes** - Google's backbone to AWS
- ✅ **Automatic validation** - Checksums verified during transfer
- ✅ **Resumable** - Survives network interruptions
- ✅ **Scheduled transfers** - Can run during off-peak hours

### How It Works

```
┌─────────────────┐         Google Transfer Service         ┌─────────────────┐
│                 │         (Google's Network Backbone)      │                 │
│  GCS Bucket(s)  │  ═══════════════════════════════════>   │   AWS S3 Bucket │
│                 │         No Egress Charges!              │                 │
└─────────────────┘                                         └─────────────────┘
        │                                                            │
        │                                                            │
        ▼                                                            ▼
  gs://motleybio/                                         s3://motleybio-data/
```

---

## Cost Analysis

### Transfer Costs (Zero! 🎉)
- **GCS → S3 via Transfer Service**: $0 (no egress charges)
- **Data transfer through Google's backbone**: $0

### AWS S3 Storage Costs (After Transfer)
Assuming 50TB of data:

| Storage Class | Monthly Cost | Best For |
|---------------|--------------|----------|
| S3 Standard | $1,150 | Frequently accessed data |
| S3 Intelligent-Tiering | $575-$1,150 | Mixed access patterns (recommended) |
| S3 Glacier Instant Retrieval | $200 | Archived data, occasional access |
| S3 Glacier Deep Archive | $50 | Long-term archives, rare access |

**Recommendation for Pipeline Data:**
- **Active analysis data**: S3 Intelligent-Tiering (~$575-$1,150/month for 50TB)
- **Archived runs**: S3 Glacier Instant Retrieval (~$200/month for 50TB)
- **Long-term backups**: S3 Glacier Deep Archive (~$50/month for 50TB)

### One-Time AWS Request Costs
- **PUT requests**: ~$0.005 per 1,000 requests
- **Example**: 1 million files = $5 in PUT requests

**Total Estimated Cost for 50TB Migration:**
- Transfer: **$0**
- Storage (first month): **$575-$1,150** (Intelligent-Tiering)
- PUT requests: **~$5-50** (depending on file count)

---

## Prerequisites

### Required Access

1. **Google Cloud Platform:**
   - Project with billing enabled
   - Owner or Editor role on GCS bucket(s)
   - `storage.transferJobs.create` permission

2. **AWS:**
   - AWS account with IAM user/role
   - S3 bucket create permission
   - Write access to destination bucket

3. **Authentication:**
   - GCP: Already authenticated (you're using the VM)
   - AWS: Access key ID and Secret access key

### Information to Gather

Before starting, collect:

- [ ] **GCS bucket names** and paths:
  ```
  gs://motleybio/
  gs://motleybio-results/
  gs://motleybio-archives/
  (list all buckets to transfer)
  ```

- [ ] **Estimated data sizes**:
  ```bash
  gsutil du -sh gs://motleybio/
  gsutil du -sh gs://motleybio-results/
  ```

- [ ] **AWS region** for S3 buckets (e.g., `us-east-1`)

- [ ] **AWS credentials** (will create in next section)

---

## AWS Setup

### Step 1: Create IAM User for Transfer Service

Create a dedicated IAM user that Google will use to write to S3:

```bash
# 1. Login to AWS Console
# 2. Go to IAM → Users → Create User
# Username: google-transfer-service

# 3. Attach inline policy (replace BUCKET_NAME with your S3 bucket name)
```

**IAM Policy JSON** (create as inline policy):

```json
{
  "Version": "2012-10-17",
  "Statement": [
    {
      "Effect": "Allow",
      "Action": [
        "s3:GetBucketLocation",
        "s3:ListBucket"
      ],
      "Resource": "arn:aws:s3:::BUCKET_NAME"
    },
    {
      "Effect": "Allow",
      "Action": [
        "s3:PutObject",
        "s3:GetObject",
        "s3:DeleteObject"
      ],
      "Resource": "arn:aws:s3:::BUCKET_NAME/*"
    }
  ]
}
```

**Replace `BUCKET_NAME`** with your actual S3 bucket name (e.g., `motleybio-data`).

### Step 2: Create Access Keys

```bash
# In AWS Console:
# 1. IAM → Users → google-transfer-service
# 2. Security credentials → Create access key
# 3. Use case: Third-party service
# 4. Download CSV or copy:
#    - Access Key ID: AKIA...
#    - Secret Access Key: wJalr...
```

**⚠️ CRITICAL:** Save these credentials securely! You'll need them in GCP Console.

### Step 3: Create S3 Destination Bucket

```bash
# Option A: Using AWS Console
# S3 → Create bucket → Choose name and region

# Option B: Using AWS CLI
aws s3 mb s3://motleybio-data --region us-east-1

# Verify bucket created
aws s3 ls s3://motleybio-data/
```

**Recommended S3 Bucket Settings:**
- **Versioning**: Enabled (allows recovery from accidental overwrites)
- **Encryption**: SSE-S3 (default encryption)
- **Block Public Access**: All enabled (keep data private)
- **Bucket Policy**: None needed (IAM user has permissions)

---

## GCP Setup

### Step 1: Enable Transfer Service API

```bash
# From GCP Console or Cloud Shell
gcloud services enable storagetransfer.googleapis.com

# Verify enabled
gcloud services list --enabled | grep storagetransfer
```

### Step 2: Grant Transfer Service Permissions

The Transfer Service needs read access to your GCS bucket:

```bash
# Get your project ID
PROJECT_ID=$(gcloud config get-value project)

# Get Transfer Service account
TRANSFER_SA=$(gcloud projects describe $PROJECT_ID \
  --format="value(projectNumber)")@storage-transfer-service.iam.gserviceaccount.com

echo "Transfer Service Account: $TRANSFER_SA"

# Grant read access to GCS bucket
gsutil iam ch serviceAccount:${TRANSFER_SA}:objectViewer gs://motleybio/
gsutil iam ch serviceAccount:${TRANSFER_SA}:legacyBucketReader gs://motleybio/

# Repeat for each bucket to transfer
```

---

## Creating the Transfer Job

### Option 1: Using GCP Console (Recommended for First Time)

This is the easiest way to set up your first transfer:

1. **Navigate to Transfer Service:**
   - Go to: https://console.cloud.google.com/transfer
   - Or: GCP Console → Cloud Storage → Transfer

2. **Create Transfer Job:**
   - Click **"Create Transfer Job"**

3. **Source Configuration:**
   - **Source type**: Google Cloud Storage
   - **Source bucket**: `motleybio` (or your bucket name)
   - **Source path (optional)**: Leave blank for entire bucket, or specify path like `results/mot26/`

4. **Destination Configuration:**
   - **Destination type**: Amazon S3
   - **AWS access key ID**: (paste from Step 2 above)
   - **AWS secret access key**: (paste from Step 2 above)
   - **S3 bucket name**: `motleybio-data`
   - **S3 bucket path (optional)**: Leave blank or specify subdirectory

5. **Transfer Options:**
   - **Overwrite objects**: "Overwrite destination if different"
   - **Delete objects**: "Don't delete objects from destination"
   - **Transfer mode**: "One-time transfer" (or schedule if needed)

6. **Settings:**
   - **Description**: "GCS to S3 Migration - motleybio bucket"
   - **Notification**: Enable Cloud Logging

7. **Schedule (for one-time transfer):**
   - **Start date/time**: Now (or schedule for off-peak)
   - **Repeat frequency**: Do not repeat

8. **Review and Create:**
   - Review all settings
   - Click **"Create"**

### Option 2: Using gcloud CLI

For automation or multiple buckets:

```bash
# Create transfer job configuration
cat > transfer-job.json <<'EOF'
{
  "description": "Transfer motleybio bucket to AWS S3",
  "status": "ENABLED",
  "projectId": "YOUR_PROJECT_ID",
  "schedule": {
    "scheduleStartDate": {
      "year": 2025,
      "month": 11,
      "day": 20
    },
    "startTimeOfDay": {
      "hours": 20,
      "minutes": 0
    }
  },
  "transferSpec": {
    "gcsDataSource": {
      "bucketName": "motleybio"
    },
    "awsS3DataSink": {
      "bucketName": "motleybio-data",
      "awsAccessKey": {
        "accessKeyId": "YOUR_AWS_ACCESS_KEY_ID",
        "secretAccessKey": "YOUR_AWS_SECRET_KEY"
      }
    },
    "objectConditions": {},
    "transferOptions": {
      "overwriteObjectsAlreadyExistingInSink": true,
      "deleteObjectsUniqueInSink": false
    }
  }
}
EOF

# Replace placeholders
sed -i 's/YOUR_PROJECT_ID/'"$PROJECT_ID"'/g' transfer-job.json
sed -i 's/YOUR_AWS_ACCESS_KEY_ID/AKIA.../g' transfer-job.json
sed -i 's/YOUR_AWS_SECRET_KEY/wJalr.../g' transfer-job.json

# Create the transfer job
gcloud transfer jobs create transfer-job.json

# List all transfer jobs
gcloud transfer jobs list
```

### Advanced Options

**Filtering specific files:**
```json
"objectConditions": {
  "includePrefixes": ["results/mot26/", "results/mot30/"],
  "excludePrefixes": ["results/temp/", "results/cache/"]
}
```

**Scheduled recurring transfers:**
```json
"schedule": {
  "scheduleStartDate": {"year": 2025, "month": 11, "day": 20},
  "scheduleEndDate": {"year": 2025, "month": 12, "day": 31},
  "startTimeOfDay": {"hours": 2, "minutes": 0},
  "repeatInterval": "86400s"  // Daily (in seconds)
}
```

---

## Transferring Specific Parts of Buckets

**One of the most powerful features:** You can transfer specific directories, exclude unwanted files, and filter by time - crucial for managing dozens of TB efficiently!

### Why Partial Transfers Matter

Pipeline outputs typically contain:
- **5-10% final results** (MultiQC, VCFs, methylation calls) ← Transfer these
- **20-30% intermediate files** (BAMs, sorted files) ← Maybe transfer
- **60-70% working files** (Nextflow work/, temp/, cache/) ← Skip these!

**Example:** A 50TB bucket with work directories can be reduced to 15TB by excluding regenerable files:
- **Transfer time**: 3 days → 18 hours (60% faster)
- **S3 storage cost**: $1,150/month → $345/month (70% savings)
- **No data loss**: Work directories are fully regenerable by pipeline

### Filtering Methods

#### **Method 1: Source Path (Single Directory)**

Transfer only a specific directory within the bucket:

**Using GCP Console:**
- Source bucket: `motleybio`
- Source path: `results/mot26/` ← Only transfers this directory

**Using CLI:**
```json
"gcsDataSource": {
  "bucketName": "motleybio",
  "path": "results/mot26/"
}
```

**Result:** Transfers `gs://motleybio/results/mot26/*` only

#### **Method 2: Include Prefixes (Multiple Specific Paths)**

Transfer multiple specific directories:

```json
"objectConditions": {
  "includePrefixes": [
    "results/mot26/",
    "results/mot30/",
    "references/",
    "important_outputs/"
  ]
}
```

**Result:** Transfers ONLY files matching these prefixes. Everything else is skipped.

**Example Use Cases:**
- Transfer specific run outputs: `"results/mot26/"`, `"results/mot30/"`
- Transfer only reference files: `"Resources/reference_genome/"`, `"Resources/GTF/"`
- Transfer only final outputs: `"results/*/multiqc/"`, `"results/*/variants/"`

#### **Method 3: Exclude Prefixes (Skip Unwanted Data)**

Transfer everything EXCEPT certain paths - **most useful for skipping work directories:**

```json
"objectConditions": {
  "excludePrefixes": [
    "work/",                    // Skip Nextflow work directories
    "temp/",                    // Skip temporary files
    ".nextflow/",               // Skip Nextflow cache
    "cache/",                   // Skip cache directories
    "logs/",                    // Skip log files
    "*/work/",                  // Skip work dirs at any level
    "*/.nextflow/",             // Skip .nextflow at any level
    "results/*/work/",          // Skip work within results
    "intermediate/",            // Skip intermediate outputs
    "scratch/"                  // Skip scratch space
  ]
}
```

**Pipeline-Specific Exclusions (Recommended):**
```json
"excludePrefixes": [
  "work/",
  ".nextflow/",
  "*/work/",
  "*/.nextflow/",
  "*.log",
  "*/trimgalore/*.fq.gz",     // Skip trimmed FASTQs (regenerable)
  "*/fastqc/*.zip",           // Skip FastQC archives (keep HTML)
  "*/star/*Aligned.out.sam",  // Skip unsorted SAM files
  "*/temp/",
  "*/cache/"
]
```

**Result:** Transfers entire bucket minus excluded paths - can reduce size by 60-80%!

#### **Method 4: Combine Include + Exclude (Most Powerful)**

Include specific areas, but exclude waste within them:

```json
"objectConditions": {
  "includePrefixes": [
    "results/mot26/",
    "results/mot30/"
  ],
  "excludePrefixes": [
    "results/mot26/work/",
    "results/mot26/temp/",
    "results/mot26/.nextflow/",
    "results/mot30/work/",
    "results/mot30/temp/",
    "results/mot30/.nextflow/"
  ]
}
```

**Result:** Transfers mot26 and mot30 outputs, but skips work/temp/cache within each - best of both worlds!

#### **Method 5: Time-Based Filtering**

Transfer only files modified within a time range:

**Recent files only (last 30 days):**
```json
"objectConditions": {
  "maxTimeElapsedSinceLastModification": "2592000s"  // 30 days in seconds
}
```

**Files modified after a specific date:**
```json
"objectConditions": {
  "lastModifiedSince": "2025-01-01T00:00:00Z",
  "lastModifiedBefore": "2025-11-20T00:00:00Z"
}
```

**Use cases:**
- Incremental transfers: Only get files created since last transfer
- Archive old data: Transfer only files older than 1 year to Glacier
- Recent results only: Skip old runs, transfer current analysis

#### **Method 6: Combine Multiple Filters**

Stack all filters for maximum control:

```json
"objectConditions": {
  "includePrefixes": [
    "results/mot26/",
    "results/mot30/"
  ],
  "excludePrefixes": [
    "results/mot26/work/",
    "results/mot30/work/"
  ],
  "lastModifiedSince": "2025-01-01T00:00:00Z"
}
```

**Result:** Transfer mot26 and mot30 (without work dirs), but only files created in 2025!

### Practical Examples for Pipeline Data

#### **Example 1: Transfer Only Final Outputs (Skip All Working Files)**

**Goal:** Transfer results, skip all intermediate/working files

```json
{
  "transferSpec": {
    "gcsDataSource": {
      "bucketName": "motleybio"
    },
    "awsS3DataSink": {
      "bucketName": "motleybio-data"
    },
    "objectConditions": {
      "excludePrefixes": [
        "work/",
        ".nextflow/",
        "cache/",
        "temp/",
        "*/work/",
        "*/.nextflow/",
        "*/temp/",
        "*/cache/",
        "intermediate/",
        "scratch/",
        "*/.command.*",
        "*.log"
      ]
    }
  }
}
```

**Potential Savings:** 50TB → 10TB (80% reduction)

#### **Example 2: Transfer Specific Run Outputs Only**

**Goal:** Transfer only QC reports, methylation data, and variants from specific runs

```json
{
  "transferSpec": {
    "gcsDataSource": {
      "bucketName": "motleybio"
    },
    "awsS3DataSink": {
      "bucketName": "motleybio-data"
    },
    "objectConditions": {
      "includePrefixes": [
        "results/mot26/multiqc/",
        "results/mot26/methylation/",
        "results/mot26/variants/",
        "results/mot26/pipeline_info/",
        "results/mot30/multiqc/",
        "results/mot30/methylation/",
        "results/mot30/variants/",
        "results/mot30/pipeline_info/"
      ]
    }
  }
}
```

**Potential Savings:** 50TB → 2TB (96% reduction) - only critical outputs!

#### **Example 3: Transfer Everything Except Work Directories**

**Goal:** Keep all outputs but skip regenerable working files

```json
{
  "transferSpec": {
    "gcsDataSource": {
      "bucketName": "motleybio"
    },
    "awsS3DataSink": {
      "bucketName": "motleybio-data"
    },
    "objectConditions": {
      "excludePrefixes": [
        "work/",
        ".nextflow/",
        "*/work/",
        "*/.nextflow/"
      ]
    }
  }
}
```

**Potential Savings:** 50TB → 20TB (60% reduction) - balanced approach

#### **Example 4: Incremental Transfer (Only New Data)**

**Goal:** Transfer only files created in the last 7 days (useful for ongoing migrations)

```json
{
  "transferSpec": {
    "gcsDataSource": {
      "bucketName": "motleybio"
    },
    "awsS3DataSink": {
      "bucketName": "motleybio-data"
    },
    "objectConditions": {
      "maxTimeElapsedSinceLastModification": "604800s",  // 7 days
      "excludePrefixes": [
        "work/",
        "*/.nextflow/"
      ]
    }
  }
}
```

**Use case:** Already transferred old data, now get new runs without re-transferring everything

#### **Example 5: Transfer References and Recent Results Only**

**Goal:** Get references + recent pipeline outputs for AWS testing

```json
{
  "transferSpec": {
    "gcsDataSource": {
      "bucketName": "motleybio"
    },
    "awsS3DataSink": {
      "bucketName": "motleybio-data"
    },
    "objectConditions": {
      "includePrefixes": [
        "Resources/",
        "results/"
      ],
      "excludePrefixes": [
        "results/*/work/",
        "results/*/.nextflow/"
      ],
      "lastModifiedSince": "2025-10-01T00:00:00Z"
    }
  }
}
```

**Potential Savings:** 50TB → 5TB - perfect for initial AWS setup!

### Using Console UI for Partial Transfers

**Step-by-step in GCP Console:**

1. **Create Transfer Job** → **Source Configuration:**
   - Source bucket: `motleybio`
   - **Source path (optional)**: Enter specific directory like `results/mot26/`

2. **Destination Configuration:**
   - Destination bucket: `motleybio-data`
   - S3 bucket path: Leave blank or specify subdirectory

3. **Click "Advanced options"** to expand filtering:

   **Include prefixes:**
   - Click "+ Add prefix"
   - Add each path to include: `results/mot26/`, `results/mot30/`, etc.
   - Only files matching these prefixes will transfer

   **Exclude prefixes:**
   - Click "+ Add prefix"
   - Add each path to exclude: `work/`, `.nextflow/`, `*/work/`, etc.
   - Files matching these will be skipped

4. **Object conditions (further down):**
   - **Created after**: Date picker for recent files
   - **Created before**: Date picker for old files

5. **Review filtering impact:**
   - Console shows estimated object count based on filters
   - Verify the count matches expectations before starting

**No JSON needed** - just fill in the form fields!

### Phased Migration Strategy for Dozens of TB

**Recommended approach for large-scale migrations:**

#### **Phase 1: References & Small Data (Immediate - Day 1)**

**Purpose:** Get AWS pipeline working quickly

```json
"objectConditions": {
  "includePrefixes": [
    "Resources/",
    "references/"
  ]
}
```

- **Size**: ~10-50GB
- **Time**: <2 hours
- **Priority**: CRITICAL (needed for AWS testing)

#### **Phase 2: Recent Important Runs (Week 1)**

**Purpose:** Transfer active analysis data

```json
"objectConditions": {
  "includePrefixes": [
    "results/mot26/",
    "results/mot30/"
  ],
  "excludePrefixes": [
    "results/mot26/work/",
    "results/mot30/work/",
    "results/*/.nextflow/"
  ]
}
```

- **Size**: 50TB → 10TB (with exclusions)
- **Time**: 12-24 hours
- **Priority**: HIGH (current projects)

#### **Phase 3: Additional Results (Week 2)**

**Purpose:** Transfer remaining valuable outputs

```json
"objectConditions": {
  "includePrefixes": [
    "results/"
  ],
  "excludePrefixes": [
    "results/mot26/",      // Already transferred
    "results/mot30/",      // Already transferred
    "results/*/work/",
    "results/*/.nextflow/",
    "results/*/temp/"
  ]
}
```

- **Size**: Variable (depends on remaining runs)
- **Time**: 1-2 days
- **Priority**: MEDIUM

#### **Phase 4: Archives (After AWS Validation)**

**Purpose:** Deep archive old data to Glacier

```json
"objectConditions": {
  "lastModifiedBefore": "2024-01-01T00:00:00Z",  // Old data
  "excludePrefixes": [
    "work/",
    "*/.nextflow/"
  ]
},
"transferOptions": {
  "s3StorageClass": "GLACIER_DEEP_ARCHIVE"  // $1/TB/month
}
```

- **Size**: Variable
- **Time**: Days (low priority)
- **Priority**: LOW
- **Cost**: Ultra-cheap ($1/TB/month vs $23/TB standard)

### Estimating Impact of Filters

**Before creating transfer job, estimate data reduction:**

```bash
# Check total bucket size
gsutil du -sh gs://motleybio/

# Check size of specific paths to include
gsutil du -sh gs://motleybio/results/mot26/
gsutil du -sh gs://motleybio/results/mot30/

# Check size of paths you're excluding (work directories)
gsutil du -sh gs://motleybio/work/
gsutil du -sh gs://motleybio/results/*/work/

# Calculate savings
# Example:
# Total: 50TB
# work/ directories: 35TB
# Expected transfer: 15TB (70% savings!)
```

**Cost/Time Comparison:**

| Scenario | Data Size | Transfer Time | S3 Cost/Month | Approach |
|----------|-----------|---------------|---------------|----------|
| No filtering | 50TB | 2-3 days | $1,150 | ❌ Wasteful |
| Exclude work/ | 15TB | 18 hours | $345 | ✅ Good |
| Only final outputs | 5TB | 6 hours | $115 | ✅ Best for testing |
| References only | 10GB | 1 hour | $0.23 | ✅ Quick start |

### Best Practices for Filtering

1. **Start small, expand later:**
   - Transfer references first (test AWS pipeline)
   - Then transfer one run completely
   - Validate before transferring everything

2. **Always exclude work directories:**
   - Work directories are 60-80% of pipeline data
   - Fully regenerable by running pipeline
   - Zero impact on reproducibility

3. **Use dry-run approach:**
   - Create transfer job
   - Check estimated object count
   - If count seems wrong, adjust filters
   - Delete and recreate job until filters are correct

4. **Document your filters:**
   - Save filter configurations
   - Note why each path is included/excluded
   - Helps with future incremental transfers

5. **Combine filters thoughtfully:**
   - Include prefixes = "only these"
   - Exclude prefixes = "except these"
   - Both together = "only these, except these parts"

6. **Plan for incrementals:**
   - After initial transfer, use time-based filters
   - Only transfer new files created since last sync
   - Saves re-transferring unchanged data

### Testing Your Filters

**Before transferring dozens of TB, test with small subset:**

```bash
# 1. Create test transfer job with filters
# 2. Use small source path to test
"gcsDataSource": {
  "bucketName": "motleybio",
  "path": "results/mot26/"  # Just one run
}

# 3. Run transfer
# 4. Validate results:
aws s3 ls s3://motleybio-data/results/mot26/ --recursive | wc -l
gsutil ls -r gs://motleybio/results/mot26/ | wc -l

# 5. If counts match and content is correct, scale up filters
```

### Common Filtering Patterns

**Skip all intermediate files:**
```json
"excludePrefixes": [
  "work/", ".nextflow/", "cache/", "temp/", "*/work/",
  "*/.nextflow/", "*/temp/", "*/cache/", "intermediate/",
  "*.log", "*/.command.*", "scratch/"
]
```

**Transfer only specific file types (using path patterns):**
```json
"includePrefixes": [
  "results/*/multiqc/multiqc_report.html",
  "results/*/variants/*.vcf.gz",
  "results/*/methylation/*.bedGraph",
  "results/*/pipeline_info/"
]
```

**Incremental daily sync:**
```json
"objectConditions": {
  "maxTimeElapsedSinceLastModification": "86400s",  // Last 24 hours
  "excludePrefixes": ["work/", "*/.nextflow/"]
},
"schedule": {
  "startTimeOfDay": {"hours": 2, "minutes": 0},
  "repeatInterval": "86400s"  // Run daily at 2 AM
}
```

---

## Monitoring the Transfer

### Using GCP Console

1. **Navigate to Transfer Jobs:**
   - https://console.cloud.google.com/transfer/jobs

2. **View Job Status:**
   - Click on your transfer job name
   - Monitor:
     - **Status**: QUEUED → IN_PROGRESS → SUCCESS
     - **Bytes transferred**: Real-time progress
     - **Objects transferred**: File count
     - **Errors**: Any failed transfers
     - **Estimated time remaining**

3. **Operation Details:**
   - Click on individual operations to see:
     - Start/end time
     - Bytes and objects transferred
     - Error summary
     - Transfer rate (MB/s)

### Using gcloud CLI

```bash
# List all transfer jobs
gcloud transfer jobs list

# Get specific job details
gcloud transfer jobs describe JOB_NAME

# Monitor active operations
gcloud transfer operations list \
  --job-names=JOB_NAME \
  --operation-statuses=in_progress

# Watch operation progress (updates every 30 seconds)
watch -n 30 'gcloud transfer operations list --job-names=JOB_NAME'
```

### Using Cloud Logging

```bash
# View transfer logs
gcloud logging read "resource.type=storage_transfer_job" \
  --limit 50 \
  --format json

# Filter for errors only
gcloud logging read "resource.type=storage_transfer_job AND severity>=ERROR" \
  --limit 50
```

### Performance Expectations

**Transfer speeds** depend on:
- Data size and object count
- Network conditions
- Geographic distance

**Typical performance:**
- **Throughput**: 1-10 Gbps (125 MB/s - 1.25 GB/s)
- **Example**: 10 TB @ 500 MB/s = ~5.5 hours
- **Example**: 50 TB @ 500 MB/s = ~28 hours
- **Example**: 100 TB @ 500 MB/s = ~56 hours (2.3 days)

**Large transfers typically run faster** due to parallelization.

---

## Validation

### Post-Transfer Verification

**1. Compare object counts:**

```bash
# GCS object count
gsutil ls -r gs://motleybio/** | wc -l

# S3 object count
aws s3 ls s3://motleybio-data/ --recursive | wc -l

# Should be equal (or very close if transfer still in progress)
```

**2. Compare total sizes:**

```bash
# GCS total size
gsutil du -sh gs://motleybio/

# S3 total size
aws s3 ls s3://motleybio-data/ --recursive --human-readable --summarize

# Should match
```

**3. Spot-check critical files:**

```bash
# Compare checksums for important files
# GCS
gsutil hash gs://motleybio/results/important_file.bam

# Download from S3 and compare
aws s3 cp s3://motleybio-data/results/important_file.bam /tmp/
md5sum /tmp/important_file.bam

# Checksums should match
```

**4. Review Transfer Service logs:**

```bash
# Check for any errors in final operation
gcloud transfer operations list \
  --job-names=JOB_NAME \
  --operation-statuses=success,failed \
  --format="table(name, status, counters.objectsCopiedToSink, counters.bytesCopiedToSink, counters.objectsFailedToDeleteFromSource)"
```

### Transfer Service Guarantees

Google Transfer Service automatically:
- ✅ Verifies checksums during transfer
- ✅ Retries failed objects (up to 3 times)
- ✅ Reports any objects that couldn't be transferred
- ✅ Maintains object metadata (when possible)

If `status=SUCCESS`, the transfer is validated by Google.

---

## Troubleshooting

### Common Issues

#### Issue 1: Permission Denied (GCS Side)

**Error:** `Permission denied when accessing GCS bucket`

**Solution:**
```bash
# Verify Transfer Service has access
gsutil iam get gs://motleybio/ | grep storage-transfer-service

# If missing, grant access
TRANSFER_SA=$(gcloud projects describe $(gcloud config get-value project) \
  --format="value(projectNumber)")@storage-transfer-service.iam.gserviceaccount.com

gsutil iam ch serviceAccount:${TRANSFER_SA}:objectViewer gs://motleybio/
gsutil iam ch serviceAccount:${TRANSFER_SA}:legacyBucketReader gs://motleybio/
```

#### Issue 2: Access Denied (AWS Side)

**Error:** `Access denied when writing to S3`

**Solution:**
```bash
# Verify IAM user policy is correct
aws iam list-user-policies --user-name google-transfer-service
aws iam get-user-policy --user-name google-transfer-service --policy-name TransferPolicy

# Test AWS credentials manually
aws s3 ls s3://motleybio-data/ \
  --profile transfer-service

# If access denied, verify:
# 1. Bucket name is correct
# 2. IAM policy includes correct bucket ARN
# 3. AWS credentials are valid and not expired
```

#### Issue 3: Transfer Stuck in QUEUED

**Error:** Job stays in QUEUED status for hours

**Solution:**
```bash
# Check if schedule is in the future
gcloud transfer jobs describe JOB_NAME --format="get(schedule)"

# Run job immediately
gcloud transfer jobs run JOB_NAME

# Or update schedule to current time
gcloud transfer jobs update JOB_NAME \
  --schedule-start-date=$(date +%Y-%m-%d) \
  --schedule-start-time=$(date +%H:%M:%S)
```

#### Issue 4: Some Objects Failed

**Error:** Operation completed with some failed objects

**Solution:**
```bash
# List failed objects
gcloud transfer operations describe OPERATION_NAME \
  --format="get(errorBreakdowns)"

# Common failures:
# - Objects deleted from source during transfer (safe to ignore)
# - Permission issues (check individual object ACLs)
# - Network timeouts (re-run job, will only copy missing objects)

# Re-run the job (automatically skips already-transferred objects)
gcloud transfer jobs run JOB_NAME
```

#### Issue 5: Slow Transfer Speed

**Problem:** Transfer is slower than expected

**Solution:**
```bash
# Check operation metrics
gcloud transfer operations describe OPERATION_NAME \
  --format="get(counters.bytesCopiedToSink, counters.bytesFoundFromSource)"

# Factors affecting speed:
# - Many small files: Slower (object overhead)
# - Few large files: Faster (better parallelization)
# - Time of day: Off-peak hours may be faster
# - S3 region: Closer to us-west (GCS primary) = faster

# For many small files, consider:
# - Using gsutil/rclone instead (better for <10GB transfers)
# - Archiving small files into tar.gz before transfer
```

---

## Best Practices

### Pre-Transfer

1. **Do a dry run with small subset:**
   ```bash
   # Transfer just one directory first
   # Set sourcePath: "results/test/" in transfer job config
   # Verify it works before transferring everything
   ```

2. **Clean up unnecessary data:**
   ```bash
   # Remove temp files, logs, caches from GCS before transfer
   gsutil rm -r gs://motleybio/temp/
   gsutil rm -r gs://motleybio/*/work/  # Nextflow work directories
   ```

3. **Document what you're transferring:**
   ```bash
   # Create inventory
   gsutil ls -lhr gs://motleybio/ > gcs_inventory_$(date +%Y%m%d).txt
   ```

4. **Choose appropriate S3 storage class:**
   - Use S3 Intelligent-Tiering for mixed access patterns
   - Set lifecycle policies for automatic archival

### During Transfer

1. **Monitor progress regularly:**
   - Check GCP Console daily
   - Set up Cloud Monitoring alerts for failures

2. **Don't modify source data:**
   - Avoid writing to GCS buckets during transfer
   - Can cause checksum mismatches or missed files

3. **Keep AWS credentials secure:**
   - Don't commit to git
   - Rotate after migration complete

### Post-Transfer

1. **Verify data integrity** (see Validation section)

2. **Update pipeline configurations:**
   ```bash
   # Update nextflow.config with S3 paths
   sed -i 's|gs://motleybio|s3://motleybio-data|g' nextflow.config
   ```

3. **Set up S3 lifecycle policies:**
   ```bash
   # Example: Archive data older than 90 days
   aws s3api put-bucket-lifecycle-configuration \
     --bucket motleybio-data \
     --lifecycle-configuration file://lifecycle-policy.json
   ```

4. **Delete GCS data after confirmed migration:**
   ```bash
   # WAIT 30+ days after transfer to ensure everything works
   # Then delete GCS buckets to stop paying for duplicate storage

   # Make final backup first!
   gsutil -m cp -r gs://motleybio/ gs://motleybio-archive-YYYYMMDD/

   # Then delete
   gsutil -m rm -r gs://motleybio/
   ```

5. **Clean up AWS IAM user:**
   ```bash
   # After transfer complete, remove transfer service IAM user
   # Or at minimum, remove/disable access keys
   aws iam delete-access-key \
     --user-name google-transfer-service \
     --access-key-id AKIA...
   ```

---

## Alternative for Specific Use Cases

### When NOT to use Google Transfer Service

**Use rclone instead if:**
- Transferring < 100 GB (rclone is faster to set up)
- Need fine-grained filtering (complex include/exclude patterns)
- Transferring specific files, not entire buckets
- Need to transform data during transfer

**Use AWS DataSync if:**
- Google Transfer Service setup is blocked by organizational policies
- Need to transfer from on-premises to S3 (not GCS→S3)
- Already have DataSync infrastructure

**Use gsutil/aws CLI if:**
- Transferring < 10 GB
- One-off manual transfer
- Need to test before automated transfer

---

## Quick Reference

### Essential Commands

```bash
# Enable Transfer Service
gcloud services enable storagetransfer.googleapis.com

# Grant permissions
TRANSFER_SA=$(gcloud projects describe $(gcloud config get-value project) \
  --format="value(projectNumber)")@storage-transfer-service.iam.gserviceaccount.com
gsutil iam ch serviceAccount:${TRANSFER_SA}:objectViewer gs://BUCKET_NAME/

# Create S3 bucket
aws s3 mb s3://DESTINATION_BUCKET --region us-east-1

# List jobs
gcloud transfer jobs list

# Monitor operation
gcloud transfer operations list --job-names=JOB_NAME

# Check for errors
gcloud transfer operations list \
  --job-names=JOB_NAME \
  --operation-statuses=failed
```

### URLs

- **GCP Transfer Console**: https://console.cloud.google.com/transfer
- **GCP Transfer Documentation**: https://cloud.google.com/storage-transfer/docs/aws-s3-transfer
- **AWS IAM Console**: https://console.aws.amazon.com/iam/
- **AWS S3 Console**: https://console.aws.amazon.com/s3/

---

## Summary

For transferring **dozens of terabytes** from GCS to S3:

✅ **Google Transfer Service saves thousands of dollars** in egress fees
✅ **Fully managed** - Google handles everything
✅ **Automatic validation** - Checksums verified
✅ **Resumable** - Survives interruptions
✅ **Fast** - Optimized Google→AWS routes

**Expected timeline for 50TB:** 1-3 days (depending on object count and network)

**Total cost:** ~$0 transfer + S3 storage (~$575-$1,150/month for Intelligent-Tiering)

---

## Need Help?

If you encounter issues not covered here:
1. Check GCP Console logs (Cloud Logging)
2. Review Transfer Service documentation: https://cloud.google.com/storage-transfer/docs
3. Contact Google Cloud Support (if you have support plan)
4. Check AWS S3 bucket policies and IAM permissions

For questions about this guide, refer to the migration checklist: `docs/gcs_to_aws_migration_checklist.md`
