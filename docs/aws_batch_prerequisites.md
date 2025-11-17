# AWS Batch Prerequisites and Setup Guide

This guide walks through everything you need to set up **before** running the methyltna pipeline on AWS Batch.

## Overview

Setting up AWS Batch involves:
1. **AWS Account Setup** (5 minutes)
2. **IAM Roles Configuration** (15 minutes)
3. **Networking Setup** (10 minutes)
4. **Batch Environment Creation** (20 minutes)
5. **S3 Bucket Setup** (5 minutes)
6. **Reference Data Upload** (varies by internet speed)
7. **Testing the Setup** (10 minutes)

**Total estimated time**: 1-2 hours (excluding reference upload time)

---

## 1. AWS Account Prerequisites

### Required
- AWS account with billing enabled
- AWS CLI installed on your local machine
- Appropriate AWS quotas for your region:
  - vCPU quota for EC2 On-Demand instances (at least 256 vCPUs)
  - vCPU quota for EC2 Spot instances (recommended: at least 256 vCPUs)

### Check Your Quotas

```bash
# Check current vCPU limits
aws service-quotas get-service-quota \
    --service-code ec2 \
    --quota-code L-1216C47A \
    --region us-east-1

# Request increase if needed (example: 512 vCPUs)
aws service-quotas request-service-quota-increase \
    --service-code ec2 \
    --quota-code L-1216C47A \
    --desired-value 512 \
    --region us-east-1
```

### Install AWS CLI (if not already installed)

```bash
# macOS
brew install awscli

# Linux
curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip"
unzip awscliv2.zip
sudo ./aws/install

# Verify installation
aws --version

# Configure with your credentials
aws configure
# Enter: Access Key ID, Secret Access Key, Region (e.g., us-east-1), Output format (json)
```

---

## 2. IAM Roles Setup

You need **three** IAM roles for AWS Batch:

### A. Batch Service Role

Allows AWS Batch to manage resources on your behalf.

```bash
# Create trust policy
cat > batch-service-trust-policy.json <<'EOF'
{
  "Version": "2012-10-17",
  "Statement": [
    {
      "Effect": "Allow",
      "Principal": {
        "Service": "batch.amazonaws.com"
      },
      "Action": "sts:AssumeRole"
    }
  ]
}
EOF

# Create the role
aws iam create-role \
    --role-name BatchServiceRole \
    --assume-role-policy-document file://batch-service-trust-policy.json

# Attach AWS managed policy
aws iam attach-role-policy \
    --role-name BatchServiceRole \
    --policy-arn arn:aws:iam::aws:policy/service-role/AWSBatchServiceRole
```

### B. ECS Task Execution Role

Allows ECS to pull Docker images and write logs.

```bash
# Create trust policy
cat > ecs-task-trust-policy.json <<'EOF'
{
  "Version": "2012-10-17",
  "Statement": [
    {
      "Effect": "Allow",
      "Principal": {
        "Service": "ecs-tasks.amazonaws.com"
      },
      "Action": "sts:AssumeRole"
    }
  ]
}
EOF

# Create the role
aws iam create-role \
    --role-name BatchECSTaskExecutionRole \
    --assume-role-policy-document file://ecs-task-trust-policy.json

# Attach AWS managed policy
aws iam attach-role-policy \
    --role-name BatchECSTaskExecutionRole \
    --policy-arn arn:aws:iam::aws:policy/service-role/AmazonECSTaskExecutionRolePolicy
```

### C. EC2 Instance Role (for Batch compute instances)

Allows EC2 instances to access S3 and ECS.

```bash
# Create trust policy
cat > ec2-trust-policy.json <<'EOF'
{
  "Version": "2012-10-17",
  "Statement": [
    {
      "Effect": "Allow",
      "Principal": {
        "Service": "ec2.amazonaws.com"
      },
      "Action": "sts:AssumeRole"
    }
  ]
}
EOF

# Create the role
aws iam create-role \
    --role-name BatchInstanceRole \
    --assume-role-policy-document file://ec2-trust-policy.json

# Attach required AWS managed policies
aws iam attach-role-policy \
    --role-name BatchInstanceRole \
    --policy-arn arn:aws:iam::aws:policy/service-role/AmazonEC2ContainerServiceforEC2Role

# Create custom S3 access policy for your pipeline bucket
cat > batch-s3-policy.json <<'EOF'
{
  "Version": "2012-10-17",
  "Statement": [
    {
      "Effect": "Allow",
      "Action": [
        "s3:GetObject",
        "s3:PutObject",
        "s3:DeleteObject",
        "s3:ListBucket"
      ],
      "Resource": [
        "arn:aws:s3:::YOUR-BUCKET-NAME/*",
        "arn:aws:s3:::YOUR-BUCKET-NAME"
      ]
    }
  ]
}
EOF

# Replace YOUR-BUCKET-NAME with your actual bucket name, then:
aws iam put-role-policy \
    --role-name BatchInstanceRole \
    --policy-name BatchS3Access \
    --policy-document file://batch-s3-policy.json

# Create instance profile
aws iam create-instance-profile \
    --instance-profile-name BatchInstanceProfile

# Add role to instance profile
aws iam add-role-to-instance-profile \
    --instance-profile-name BatchInstanceProfile \
    --role-name BatchInstanceRole
```

### D. Spot Fleet Role (for Spot instances - recommended for cost savings)

```bash
# Create trust policy
cat > spot-fleet-trust-policy.json <<'EOF'
{
  "Version": "2012-10-17",
  "Statement": [
    {
      "Effect": "Allow",
      "Principal": {
        "Service": "spotfleet.amazonaws.com"
      },
      "Action": "sts:AssumeRole"
    }
  ]
}
EOF

# Create the role
aws iam create-role \
    --role-name BatchSpotFleetRole \
    --assume-role-policy-document file://spot-fleet-trust-policy.json

# Attach AWS managed policy
aws iam attach-role-policy \
    --role-name BatchSpotFleetRole \
    --policy-arn arn:aws:iam::aws:policy/service-role/AmazonEC2SpotFleetTaggingRole
```

---

## 3. Networking Setup

AWS Batch needs a VPC, subnets, and security groups.

### Option A: Use Default VPC (Simplest)

```bash
# Get your default VPC ID
DEFAULT_VPC=$(aws ec2 describe-vpcs \
    --filters "Name=isDefault,Values=true" \
    --query "Vpcs[0].VpcId" \
    --output text)

echo "Default VPC: $DEFAULT_VPC"

# Get default subnets
aws ec2 describe-subnets \
    --filters "Name=vpc-id,Values=$DEFAULT_VPC" \
    --query "Subnets[*].SubnetId" \
    --output text

# Create security group
SECURITY_GROUP=$(aws ec2 create-security-group \
    --group-name methyltna-batch-sg \
    --description "Security group for methylTNA Batch jobs" \
    --vpc-id $DEFAULT_VPC \
    --query "GroupId" \
    --output text)

echo "Security Group: $SECURITY_GROUP"

# Allow outbound internet access (needed to pull Docker images)
aws ec2 authorize-security-group-egress \
    --group-id $SECURITY_GROUP \
    --ip-permissions IpProtocol=-1,IpRanges='[{CidrIp=0.0.0.0/0}]'
```

### Option B: Create New VPC (More Control)

<details>
<summary>Click to expand VPC creation steps</summary>

```bash
# Create VPC
VPC_ID=$(aws ec2 create-vpc \
    --cidr-block 10.0.0.0/16 \
    --query Vpc.VpcId \
    --output text)

aws ec2 create-tags \
    --resources $VPC_ID \
    --tags Key=Name,Value=methyltna-batch-vpc

# Enable DNS
aws ec2 modify-vpc-attribute \
    --vpc-id $VPC_ID \
    --enable-dns-support

aws ec2 modify-vpc-attribute \
    --vpc-id $VPC_ID \
    --enable-dns-hostnames

# Create Internet Gateway
IGW_ID=$(aws ec2 create-internet-gateway \
    --query InternetGateway.InternetGatewayId \
    --output text)

aws ec2 attach-internet-gateway \
    --vpc-id $VPC_ID \
    --internet-gateway-id $IGW_ID

# Create public subnet
SUBNET_ID=$(aws ec2 create-subnet \
    --vpc-id $VPC_ID \
    --cidr-block 10.0.1.0/24 \
    --availability-zone us-east-1a \
    --query Subnet.SubnetId \
    --output text)

# Create route table
ROUTE_TABLE_ID=$(aws ec2 create-route-table \
    --vpc-id $VPC_ID \
    --query RouteTable.RouteTableId \
    --output text)

aws ec2 create-route \
    --route-table-id $ROUTE_TABLE_ID \
    --destination-cidr-block 0.0.0.0/0 \
    --gateway-id $IGW_ID

aws ec2 associate-route-table \
    --subnet-id $SUBNET_ID \
    --route-table-id $ROUTE_TABLE_ID

# Create security group (same as Option A)
SECURITY_GROUP=$(aws ec2 create-security-group \
    --group-name methyltna-batch-sg \
    --description "Security group for methylTNA Batch jobs" \
    --vpc-id $VPC_ID \
    --query "GroupId" \
    --output text)
```
</details>

---

## 4. Create AWS Batch Resources

### A. Compute Environment

```bash
# Get your subnet IDs (use default VPC or the ones you created)
SUBNET_IDS=$(aws ec2 describe-subnets \
    --filters "Name=vpc-id,Values=$DEFAULT_VPC" \
    --query "Subnets[*].SubnetId" \
    --output text | tr '\t' ',')

# Get role ARNs
BATCH_SERVICE_ROLE=$(aws iam get-role --role-name BatchServiceRole --query Role.Arn --output text)
INSTANCE_PROFILE=$(aws iam get-instance-profile --instance-profile-name BatchInstanceProfile --query InstanceProfile.Arn --output text)
SPOT_FLEET_ROLE=$(aws iam get-role --role-name BatchSpotFleetRole --query Role.Arn --output text)

# Create compute environment (SPOT instances for cost savings)
cat > compute-environment.json <<EOF
{
  "computeEnvironmentName": "methyltna-spot-compute",
  "type": "MANAGED",
  "state": "ENABLED",
  "computeResources": {
    "type": "SPOT",
    "allocationStrategy": "SPOT_CAPACITY_OPTIMIZED",
    "minvCpus": 0,
    "maxvCpus": 256,
    "desiredvCpus": 0,
    "instanceTypes": [
      "m5.2xlarge",
      "m5.4xlarge",
      "r5.2xlarge",
      "r5.4xlarge",
      "c5.4xlarge"
    ],
    "subnets": $(echo $SUBNET_IDS | jq -R 'split(",")'),
    "securityGroupIds": ["$SECURITY_GROUP"],
    "instanceRole": "$INSTANCE_PROFILE",
    "spotIamFleetRole": "$SPOT_FLEET_ROLE",
    "bidPercentage": 100,
    "tags": {
      "Name": "methyltna-batch-instance",
      "Project": "methylTNA"
    }
  },
  "serviceRole": "$BATCH_SERVICE_ROLE"
}
EOF

aws batch create-compute-environment --cli-input-json file://compute-environment.json
```

**Wait for compute environment to be VALID:**

```bash
# Check status (should say VALID when ready)
aws batch describe-compute-environments \
    --compute-environments methyltna-spot-compute \
    --query "computeEnvironments[0].status"
```

### B. Job Queue

```bash
# Get compute environment ARN
COMPUTE_ENV_ARN=$(aws batch describe-compute-environments \
    --compute-environments methyltna-spot-compute \
    --query "computeEnvironments[0].computeEnvironmentArn" \
    --output text)

# Create job queue
cat > job-queue.json <<EOF
{
  "jobQueueName": "methyltna-queue",
  "state": "ENABLED",
  "priority": 1,
  "computeEnvironmentOrder": [
    {
      "order": 1,
      "computeEnvironment": "$COMPUTE_ENV_ARN"
    }
  ]
}
EOF

aws batch create-job-queue --cli-input-json file://job-queue.json
```

---

## 5. S3 Bucket Setup

```bash
# Create bucket (replace UNIQUE-BUCKET-NAME with your desired name)
BUCKET_NAME="methyltna-pipeline-YOUR-NAME"

aws s3 mb s3://$BUCKET_NAME --region us-east-1

# Create directory structure
aws s3api put-object --bucket $BUCKET_NAME --key work/
aws s3api put-object --bucket $BUCKET_NAME --key results/
aws s3api put-object --bucket $BUCKET_NAME --key references/
aws s3api put-object --bucket $BUCKET_NAME --key input/

# Optional: Set lifecycle policy to auto-delete old work files
cat > lifecycle-policy.json <<'EOF'
{
  "Rules": [
    {
      "Id": "DeleteOldWorkFiles",
      "Status": "Enabled",
      "Prefix": "work/",
      "Expiration": {
        "Days": 30
      }
    }
  ]
}
EOF

aws s3api put-bucket-lifecycle-configuration \
    --bucket $BUCKET_NAME \
    --lifecycle-configuration file://lifecycle-policy.json

# Enable versioning (optional but recommended)
aws s3api put-bucket-versioning \
    --bucket $BUCKET_NAME \
    --versioning-configuration Status=Enabled
```

---

## 6. Upload Reference Files to S3

### Upload Genome and Annotations

```bash
# Upload genome FASTA
aws s3 cp /path/to/genome.fa s3://$BUCKET_NAME/references/genome.fa

# Upload GTF annotations
aws s3 cp /path/to/annotations.gtf s3://$BUCKET_NAME/references/annotations.gtf
```

### Option A: Upload Pre-built Indexes (Recommended - saves 2 hours per run)

```bash
# Upload STAR index directory
aws s3 sync /path/to/star_index/ s3://$BUCKET_NAME/references/star_index/

# Upload Biscuit index directory
aws s3 sync /path/to/biscuit_index/ s3://$BUCKET_NAME/references/biscuit_index/
```

### Option B: Let Pipeline Build Indexes (will take ~2 hours on first run)

If you skip uploading indexes, the pipeline will build them on first run and you can save them for reuse:

```bash
# After first run completes, download and save indexes
aws s3 sync s3://$BUCKET_NAME/work/star_genome_generate/ ./star_index/
aws s3 sync s3://$BUCKET_NAME/work/biscuit_index/ ./biscuit_index/

# Upload for future runs
aws s3 sync ./star_index/ s3://$BUCKET_NAME/references/star_index/
aws s3 sync ./biscuit_index/ s3://$BUCKET_NAME/references/biscuit_index/
```

---

## 7. Test Your Setup

### Create Test Samplesheet

```bash
cat > test-samplesheet.csv <<'EOF'
sample,fastq_1,fastq_2
test_sample,s3://nf-core-awsmegatests/methylseq/input/SRR389222_sub1.fastq.gz,s3://nf-core-awsmegatests/methylseq/input/SRR389222_sub2.fastq.gz
EOF

# Upload to S3
aws s3 cp test-samplesheet.csv s3://$BUCKET_NAME/input/test-samplesheet.csv
```

### Run Test Pipeline

```bash
nextflow run Motleybio-organization/methyltna \
    -profile test,awsbatch \
    --aws_queue methyltna-queue \
    --aws_region us-east-1 \
    --outdir s3://$BUCKET_NAME/test-results \
    -w s3://$BUCKET_NAME/work
```

### Monitor Test Run

```bash
# Watch job queue
watch -n 10 'aws batch list-jobs \
    --job-queue methyltna-queue \
    --job-status RUNNING \
    --query "jobSummaryList[*].[jobName,status]" \
    --output table'

# Check compute environment scaling
aws batch describe-compute-environments \
    --compute-environments methyltna-spot-compute \
    --query "computeEnvironments[0].computeResources.[minvCpus,desiredvCpus,maxvCpus]"
```

---

## 8. Cost Optimization Checklist

Before running production workloads:

- [ ] **Use Spot instances** (compute environment set to SPOT) - saves 70-90%
- [ ] **Set maxvCpus appropriately** - don't over-provision
- [ ] **Enable S3 lifecycle policies** - auto-delete old work directories after 30 days
- [ ] **Use S3 Intelligent Tiering** - already configured in pipeline
- [ ] **Upload pre-built indexes** - saves ~2 hours of compute per run
- [ ] **Enable CloudWatch alarms** - get notified of unexpected costs
- [ ] **Tag all resources** - for cost tracking and reporting

### Estimated Costs (US-East-1, Spot instances)

For a typical 12-sample TNA-seq run:
- **Compute**: ~$20-30 (with Spot instances)
- **S3 Storage**: ~$5-10/month for work files, $1-2/month for results
- **Data Transfer**: Minimal if using S3 for input/output

**Total**: ~$25-40 per 12-sample run (vs $200-300 with On-Demand instances)

---

## Troubleshooting

### "Unable to assume role"
```bash
# Check role trust policies are correct
aws iam get-role --role-name BatchServiceRole --query Role.AssumeRolePolicyDocument
```

### "Compute environment stuck in INVALID"
```bash
# Check detailed status
aws batch describe-compute-environments \
    --compute-environments methyltna-spot-compute \
    --query "computeEnvironments[0].statusReason"

# Common issues:
# - Invalid subnet IDs
# - Security group in wrong VPC
# - Missing IAM permissions
```

### "Job fails immediately"
```bash
# Check CloudWatch logs
aws logs tail /aws/batch/job --follow

# Common issues:
# - S3 permissions
# - Invalid container image
# - Missing input files
```

---

## Cleanup (when you're done testing)

```bash
# Delete job queue
aws batch update-job-queue \
    --job-queue methyltna-queue \
    --state DISABLED

aws batch delete-job-queue \
    --job-queue methyltna-queue

# Delete compute environment
aws batch update-compute-environment \
    --compute-environment methyltna-spot-compute \
    --state DISABLED

aws batch delete-compute-environment \
    --compute-environment methyltna-spot-compute

# Delete S3 bucket (if you want)
aws s3 rb s3://$BUCKET_NAME --force

# Delete security group, VPC, etc. (if created for this purpose)
```

---

## Next Steps

Once setup is complete:
1. Read [AWS Batch Setup Guide](aws_batch_setup.md) for production usage
2. Review [main README](../README.md) for pipeline parameters
3. Start processing your samples!

## Support

- AWS Batch Documentation: https://docs.aws.amazon.com/batch/
- Nextflow AWS Batch: https://www.nextflow.io/docs/latest/awscloud.html
- Pipeline Issues: https://github.com/Motleybio-organization/methyltna/issues
