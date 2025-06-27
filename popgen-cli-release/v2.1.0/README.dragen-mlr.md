# DRAGEN PopGen Machine Learning Recalibration (MLR) Workflow

This is a brief README on how to install, configure, and launch DRAGEN MLR analyses. More detailed documentation will be provided later.  
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
tar xvzf popgen-dragen-mlr-demo.tar.gz
mv popgen-dragen-mlr-demo/* .
rm -rf popgen-dragen-mlr-demo
```

---

## Step 3: Configure the DRAGEN MLR Workflow

```bash
./1.config.sh
```

---

## Step 4: Launch MLR Analyses (DRAGEN v3.7.8, GIAB Sample gVCFs)

```bash
./2.submit.sh
```

---

## Step 5: Monitor Running Jobs

```bash
watch -d ./poll.sh
```
