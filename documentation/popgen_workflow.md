# Population Genotyping with ICA IGG

[Illumina Dragen documentation](https://help.dragen.illumina.com/product-guides/dragen-v4.5/dragen-dna-pipeline/iterative-gvcf-genotyper)

This repo contains nextflow workflows for running the iterative genotyper. These automate the manual process as described in the demo project. Input files and results are removed from ICA to minimise costs.

# Prerequisites

These steps should be carried out in a directory on the cluster that is shared with the task runner.

1. Clone the [ica-elwazi](https://github.com/SBIMB/ica-elwazi/) repository and change directory to `nextflow_workflows/popgen`.
2. Install the `popgen-cli` contained in the repository. Or get the latest from [ICA](https://help.dragen.illumina.com/dragen-v4.5/reference/release-notes-readme/dragen-crn-v4.5.4#iterative-gvcf-genotyper-igg).

```bash
cd ica-igg-pipeline

# Using Conda for support in nextflow.
conda create -p ./.conda/popgen python
conda activate .conda/popgen 
pip install Downloads/popgen_cli-2.2.3-py3-none-any.whl

# On the cluster
/opt/exp_soft/miniconda3/bin/conda create -p ./.conda/popgen python
source /opt/exp_soft/miniconda3/bin/activate .conda/popgen
module load python/3.8
pip install Downloads/popgen_cli-2.2.3-py3-none-any.whl
```

3. Install the [ICA CLI](https://help.ica.illumina.com/command-line-interface/cli-releasehistory#releases).
4. [Authenticate](https://help.ica.illumina.com/command-line-interface/cli-authentication) the CLI with your ICA [API key](https://help.ica.illumina.com/get-started/gs-getstarted#api-keys).
5. Create a project in [ICA](https://ica.illumina.com/ica/projects/) for the analyses. The project name must end with “-jobs” and only contain lower case characters, digits and dashes “-”. 
6. [Link](https://help.ica.illumina.com/home/h-bundles#linking-an-existing-bundle-to-a-project) the Popgen CLI Bundle to your project in the ICA GUI. This can be done while creating the project. 
7. Use the Popgen CLI to create and configure your local directory. You can [retrieve](https://app.notion.com/p/Command-Line-35e80cd41dae8024be69c8fdf9117907?pvs=21) your API key with `icav2 config get`. 

```bash
# venv: source .venv/bin/activate
source /opt/exp_soft/miniconda3/bin/activate .conda/popgen
popgen-cli dragen-igg config
```

8. Optional, for convenience in subsequent steps. Create a file called `project-config.json` to avoid having to pass the config file as an argument to the workflow. If you do this remember to replace it when config changes.

```bash
cp $(ls project-data/secret/project-config.*.json | tail -n1) project-data/secret/project-config.json 
```

# Run the workflow

Run screen or tmux.

Activate the conda environment.

The default path for sequence files is `./test_data/version-n` you can symlink to the data.

## Initial batch

```bash
# source /opt/exp_soft/miniconda3/bin/activate .conda/popgen
nextflow run -profile slurm workflows/initial_analyses.nf
```

## Subsequent batches

The subsequent analyses pipeline takes the results directory as input and re-uploads the census files for iterative analysis.

```bash
nextflow run -profile slurm workflows/subsequent_analyses.nf 
```

