# GCS to S3 Migration Report

**Migration Date:** November 22, 2025
**Source:** Google Cloud Storage (GCS) - `gcs:motleybio`
**Destination:** Amazon S3 - `s3:motleybio`
**Server:** rclone-s3
**Region:** us-east-2

---

## Executive Summary

Successfully migrated approximately **3+ TB** of data from Google Cloud Storage to Amazon S3, comprising genomics data, reference files, software packages, and datasets. All transfers completed successfully with checksum verification enabled.

---

## Migration Details

### 1. Laboratory Assay Development Data
**Path:** `Laboratory/ASSAY_DEVELOPMENT/`
**Total Duration:** ~12.5 hours
**Start Time:** November 22, 2025 00:58:38

| Directory | Size | Files | Duration | Avg Speed | Status |
|-----------|------|-------|----------|-----------|--------|
| concat_bams | 192.7 GB | 9 | 1h 0m | 54.6 MiB/s | ✓ Complete |
| run_motley21 | 332.1 GB | 1,104 | 1h 20m | 62.7 MiB/s | ✓ Complete |
| run_motley22 | 11.9 GB | 288 | 1m 40s | 113.4 MiB/s | ✓ Complete |
| run_motley23 | 49.4 GB | 440 | 10m 57s | 67.9 MiB/s | ✓ Complete |
| run_motley26 | 1.74 TB | 1,398 | 7h 19m | 56.6 MiB/s | ✓ Complete |
| run_motley28 | 196.2 GB | 1,851 | 47m 47s | 62.8 MiB/s | ✓ Complete |
| run_motley29 | 221.1 GB | 2,321 | 53m 39s | 69.8 MiB/s | ✓ Complete |
| run_motley30 | 238.4 GB | 1,956 | 57m 40s | 64.5 MiB/s | ✓ Complete |
| **Subtotal** | **~2.93 TB** | **9,367 files** | **~12h 32m** | **~64.7 MiB/s avg** | ✓ Complete |

**Transfer Command:**
```bash
for dir in concat_bams run_motley21 run_motley22 run_motley23 run_motley26 run_motley28 run_motley29 run_motley30; do
  echo "=== Transferring $dir ==="
  rclone copy \
    gcs:motleybio/Laboratory/ASSAY_DEVELOPMENT/$dir/ \
    s3:motleybio/Laboratory/ASSAY_DEVELOPMENT/$dir/ \
    --progress --stats 30s --transfers 16 --checkers 32 --checksum \
    --log-file ~/rclone-$dir-$(date +%Y%m%d-%H%M%S).log
done
```

---

### 2. Reference Files and Indexes
**Completion Time:** November 22, 2025 01:00-01:01

| Item | Source Path | Destination Path | Size | Files | Duration | Status |
|------|------------|------------------|------|-------|----------|--------|
| Resources | `Resources/` | `Resources/` | 11.6 MB | 259 | 3.2s | ✓ Complete |
| Methylation References | `Workspaces/Marcus/methyltna_references/` | `Resources/methyltna_references/` | 55.7 GB | 52 | 30m 51s | ✓ Complete |
| **Subtotal** | | | **55.7 GB** | **311 files** | **~31m** | ✓ Complete |

---

### 3. Software Packages
**Path:** `NGS_Software/`
**Completion Time:** November 22, 2025 19:47

| Item | Size | Files | Duration | Avg Speed | Status |
|------|------|-------|----------|-----------|--------|
| NGS_Software | 7.4 GB | 45,144 | 7m 6s | 394.9 KiB/s | ✓ Complete |

**Transfer Command:**
```bash
rclone copy \
  gcs:motleybio/NGS_Software/ \
  s3:motleybio/NGS_Software/ \
  --progress --stats 30s --transfers 16 --checkers 32 --checksum \
  --log-file ~/rclone-ngs-software-$(date +%Y%m%d-%H%M%S).log
```

---

### 4. Machine Learning Datasets
**Path:** `Datasets/`
**Completion Time:** November 22, 2025 20:11

| Dataset | Size | Files | Status |
|---------|------|-------|--------|
| FOUNDATION_MODEL_DATASETS | 592.7 GB | 47,220 | ✓ Complete |
| SEED_DATA_PACKAGE | 16.0 GB | 161 | ✓ Complete |
| **Subtotal** | **608.7 GB** | **47,381 files** | ✓ Complete |

---

## Assessed but Not Yet Transferred

The following datasets were assessed for size and content but have not yet been migrated:

| Workspace/Dataset | Size | Files | Notes |
|-------------------|------|-------|-------|
| Workspaces/Sunil/ | 1.4 TB | 10,940 | Pending transfer decision |

---

## Migration Statistics Summary

| Category | Total Size | Total Files | Status |
|----------|-----------|-------------|--------|
| Laboratory Data | ~2.93 TB | 9,367 | ✓ Complete |
| References/Resources | 55.7 GB | 311 | ✓ Complete |
| Software | 7.4 GB | 45,144 | ✓ Complete |
| ML Datasets | 608.7 GB | 47,381 | ✓ Complete |
| **Grand Total** | **~3.58 TB** | **102,203 files** | ✓ Complete |

---

## Technical Configuration

### rclone Parameters Used
- `--progress`: Display real-time progress
- `--stats 30s`: Update statistics every 30 seconds
- `--transfers 16`: Perform 16 parallel file transfers
- `--checkers 32`: Use 32 parallel checksum verification threads
- `--checksum`: Verify transfers using checksums (ensures data integrity)
- `--log-file`: Individual log file for each transfer operation

### S3 Configuration
- **Region:** us-east-2 (automatically detected and switched from us-east-1)
- **Bucket:** motleybio

---

## Log Files

All transfer operations generated detailed log files stored in `/home/marcus/`:

```
rclone-concat_bams-20251122-005838.log
rclone-run_motley21-20251122-015839.log
rclone-run_motley22-20251122-031912.log
rclone-run_motley23-20251122-032052.log
rclone-run_motley26-20251122-033149.log
rclone-run_motley28-20251122-105134.log
rclone-run_motley29-20251122-113921.log
rclone-run_motley30-20251122-123300.log
rclone-resources-20251122-010048.log
rclone-indexes-20251122-010126.log
rclone-ngs-software-20251122-194756.log
rclone-foundation-20251122-201118.log
rclone-seed-data-20251122-201117.log
```

---

## Notes and Observations

1. **Successful Completion:** All transfers completed successfully with no errors reported in log files
2. **Data Integrity:** Checksum verification was enabled for all transfers, ensuring data integrity
3. **Performance:** Average transfer speed of ~64.7 MiB/s for large genomics datasets
4. **Parallel Processing:** Used 16 parallel transfers and 32 checkers for optimal performance
5. **Region Switching:** S3 automatically detected and switched to us-east-2 region
6. **Largest Transfer:** run_motley26 (1.74 TB) took the longest at 7h 19m

---

## Next Steps / Pending Items

- [ ] Review Sunil workspace (1.4 TB) for potential migration
- [ ] Verify data accessibility in S3
- [ ] Update application configurations to use S3 paths
- [ ] Consider archiving or deleting GCS data after verification period
- [ ] Document any application-specific migration requirements

---

**Report Generated:** November 23, 2025
**Prepared By:** Migration automation on rclone-s3 server
