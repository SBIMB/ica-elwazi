# DRAGEN PopGen Iterative gVCF Genotyper (IGG) Workflow

This is a brief README on how to install, configure, and launch DRAGEN IGG analyses. More detailed documentation will be provided later.  
For any questions, please contact Zhuoyi Huang at <zhuang@illumina.com>.

---

## Step 1: Install PopGen CLI

```bash
python3 -m venv venv
source venv/bin/activate
pip install popgen_cli-2.0.0-py3-none-any.whl
popgen-cli --version
```

---

## Step 2: Unpack the Demo Package

```bash
tar xvzf popgen-dragen-igg-demo.tar.gz
mv popgen-dragen-igg-demo/* .
rm -rf popgen-dragen-igg-demo
```

---

## Step 3: Configure the DRAGEN IGG Workflow

```bash
./1.config.sh
```

---

## Step 4: Upload Demo Project Data to ICA

```bash
./2.upload.sh
```

---

## Step 5: Launch IGG Analyses (DRAGEN v4.3.13f, GIAB Sample gVCFs)

For demo purposes, analyses are restricted to chromosome 20 (`chrom-20`), which includes:
- `shard-88`
- `shard-89`

Seven GIAB samples are split into two batches:
- `batch-1`: 4 samples  
- `batch-2`: 3 samples

To simulate N+1 aggregation, two analysis versions are defined:
- **Version 1** includes only `batch-1`
- **Version 2** includes both `batch-1` and `batch-2`  
  *(no need to re-run step1 for `batch-1` in version 2)*

Follow the steps below in order for each version, ensuring all previous jobs are completed before moving forward.

---

### Version 1: Joint Analysis (batch-1 only)

```bash
version_id=1
batch_ids=1
shard_ids=88-89
chrom_ids=20

./3.submit.sh step1 $batch_ids $shard_ids
./3.submit.sh step2 $version_id $shard_ids
./3.submit.sh count_subshards $version_id $shard_ids
./3.submit.sh step3 $version_id $shard_ids
./3.submit.sh step4 $version_id $chrom_ids
```

---

### Version 2: Joint Analysis (batch-1 + batch-2)

```bash
version_id=2
batch_ids=2
shard_ids=88-89
chrom_ids=20

./3.submit.sh step1 $batch_ids $shard_ids
./3.submit.sh step2 $version_id $shard_ids
./3.submit.sh count_subshards $version_id $shard_ids
./3.submit.sh step3 $version_id $shard_ids
./3.submit.sh step4 $version_id $chrom_ids
```

---

## Step 6: Monitor Running Jobs

```bash
watch -d ./poll.sh
```
