# AWS Batch Setup Guide

This guide explains how to run the methyltna pipeline on AWS Batch.

## Prerequisites

### 1. AWS Infrastructure

You need the following AWS resources configured:

- **AWS Batch Compute Environment** with appropriate instance types
- **AWS Batch Job Queue** connected to the compute environment
- **S3 Bucket** for storing work files and results
- **IAM Roles**:
  - Batch job execution role
  - EC2 instance role with S3 access
  - Nextflow execution role (if running from EC2/ECS)

### 2. Local Setup

- **AWS CLI** installed and configured (`aws configure`)
- **Nextflow** version ≥24.10.5
- AWS credentials with permissions to:
  - Submit Batch jobs
  - Read/write to S3 bucket
  - Describe Batch resources

## Quick Start

### Basic Usage

```bash
nextflow run Motleybio-organization/methyltna \
    -profile awsbatch \
    --aws_queue my-batch-queue \
    --input s3://my-bucket/samplesheet.csv \
    --outdir s3://my-bucket/results \
    -w s3://my-bucket/work \
    --genome_fasta s3://my-bucket/references/genome.fa \
    --annotation_gtf s3://my-bucket/references/annotations.gtf
```

### With Pre-built Indexes (Recommended)

```bash
nextflow run Motleybio-organization/methyltna \
    -profile awsbatch \
    --aws_queue my-batch-queue \
    --aws_region us-east-1 \
    --input s3://my-bucket/samplesheet.csv \
    --outdir s3://my-bucket/results \
    -w s3://my-bucket/work \
    --genome_fasta s3://my-bucket/references/genome.fa \
    --annotation_gtf s3://my-bucket/references/annotations.gtf \
    --star_index s3://my-bucket/references/star_index/ \
    --biscuit_index s3://my-bucket/references/biscuit_index/
```

### Resume Failed Run

```bash
nextflow run Motleybio-organization/methyltna \
    -profile awsbatch \
    --aws_queue my-batch-queue \
    --input s3://my-bucket/samplesheet.csv \
    --outdir s3://my-bucket/results \
    -w s3://my-bucket/work \
    -resume
```

## Important Configuration

### Required Parameters

| Parameter | Description | Example |
|-----------|-------------|---------|
| `--aws_queue` | AWS Batch job queue name | `tna-seq-queue` |
| `-w` | S3 work directory | `s3://bucket/work` |
| `--input` | Samplesheet (can be S3 or local) | `s3://bucket/samples.csv` |
| `--outdir` | Results directory (can be S3 or local) | `s3://bucket/results` |

### Optional Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--aws_region` | AWS region | `us-east-1` |

### Resource Optimization

The pipeline is configured with sensible defaults for AWS Batch:

- **STAR alignment**: 64 GB memory
- **Biscuit alignment**: 32 GB memory
- **RSeQC gene body coverage**: 16 GB memory (with 2M subsampling)
- **Standard processes**: 8 GB memory

You can override these in a custom config file if needed.

## Cost Optimization Tips

1. **Use Pre-built Indexes**
   - Upload STAR and Biscuit indexes to S3
   - Specify with `--star_index` and `--biscuit_index`
   - Saves ~2 hours of compute per run

2. **Use Spot Instances**
   - Configure your Batch compute environment with EC2 Spot
   - Can save 70-90% on compute costs
   - Pipeline is fault-tolerant with retry logic

3. **Enable RSeQC Subsampling** (Default)
   - Automatically enabled with 2M read subsampling
   - Reduces RSeQC time by 8-10x
   - Disable with `--skip_rseqc_subsampling` if full-depth QC needed

4. **S3 Storage Class**
   - Results use INTELLIGENT_TIERING by default
   - Work directory can be deleted after successful completion
   - Or set lifecycle policy to delete after 30 days

## Monitoring

### View Running Jobs

```bash
aws batch list-jobs \
    --job-queue my-batch-queue \
    --job-status RUNNING
```

### Check Job Logs

```bash
aws logs tail /aws/batch/job \
    --follow \
    --log-stream-name <job-id>
```

### Nextflow Tower (Optional)

For better monitoring, use Nextflow Tower:

```bash
export TOWER_ACCESS_TOKEN="your-token"

nextflow run Motleybio-organization/methyltna \
    -profile awsbatch \
    --aws_queue my-batch-queue \
    -with-tower
```

## Troubleshooting

### Common Issues

**Error: "AWS Batch queue not specified"**
```bash
# Solution: Add --aws_queue parameter
--aws_queue my-batch-queue
```

**Error: "S3 work directory required"**
```bash
# Solution: Specify S3 work directory
-w s3://my-bucket/work
```

**Error: "Container pull timeout"**
```bash
# Solution: Increase timeout in your Batch compute environment settings
# Or use pre-pulled containers in ECR
```

**High Costs**
```bash
# Solutions:
# 1. Use Spot instances in compute environment
# 2. Enable auto-scaling to terminate idle instances
# 3. Set --keep_intermediates false (default)
# 4. Use S3 lifecycle policies to delete old work directories
```

## Example Batch Compute Environment Setup

Here's a Terraform example for reference:

```hcl
resource "aws_batch_compute_environment" "methyltna" {
  compute_environment_name = "methyltna-compute"
  type                     = "MANAGED"

  compute_resources {
    type                = "SPOT"
    allocation_strategy = "SPOT_CAPACITY_OPTIMIZED"

    instance_types = [
      "m5.2xlarge",   # 8 vCPU, 32 GB - general purpose
      "m5.4xlarge",   # 16 vCPU, 64 GB - STAR alignment
      "r5.2xlarge",   # 8 vCPU, 64 GB - memory intensive
      "r5.4xlarge",   # 16 vCPU, 128 GB - high memory
    ]

    min_vcpus     = 0
    desired_vcpus = 0
    max_vcpus     = 256

    subnets         = var.subnet_ids
    security_groups = [var.security_group_id]

    instance_role = aws_iam_instance_profile.batch_instance.arn

    spot_iam_fleet_role = aws_iam_role.spot_fleet.arn
    bid_percentage      = 100
  }

  service_role = aws_iam_role.batch_service.arn
}

resource "aws_batch_job_queue" "methyltna" {
  name     = "methyltna-queue"
  state    = "ENABLED"
  priority = 1

  compute_environments = [
    aws_batch_compute_environment.methyltna.arn
  ]
}
```

## Support

For issues specific to AWS Batch:
- Check CloudWatch logs for Batch jobs
- Review Nextflow trace reports in results directory
- Ensure IAM permissions are correct

For pipeline issues:
- See main README.md
- Check GitHub issues: https://github.com/Motleybio-organization/methyltna/issues
