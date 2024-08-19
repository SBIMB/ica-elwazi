## Run Nextflow Pipelines (Analyses)   
We'd like to perform an analysis on a file using the CLI. For convenience and for testing purposes, we will use small files. We can download small files of various formats from [here](https://ftp.ncbi.nih.gov/genomes/HUMAN_MICROBIOM/Bacteria/). Since the Nextflow pipeline that exists in our project performs an analysis on FASTA files, we'll use the FASTA format for our tests.   

We can use the UI to kick off an analysis on an uploaded file. The analysis from the example in the tutorial takes about 30 minutes to complete. When completed, there is output data in the ICA storage. We can get a list of analyses belonging to a project with the following command:
```bash
icav2 projectanalyses list --project-id <projectId>
```

We can filter by `userReference` to get the `id` or `status` of the desired analysis. The CLI can be used to get the status of the analysis:
```bash
icav2 projectanalyses get <analysisId>
```
The JSON response contains a lot of information, but we are interested in the status of the analysis, i.e.
```json
{
	"analysisPriority": "MEDIUM",
	"analysisStorage": {
		"description": "1.2TB",
		"id": "6e1b6c8f-f913-48b2-9bd0-7fc13eda0fd0",
		"name": "Small",
		"ownerId": "8ec463f6-1acb-341b-b321-043c39d8716a",
		"tenantId": "f91bb1a0-c55f-4bce-8014-b2e60c0ec7d3",
		"tenantName": "ica-cp-admin",
		"timeCreated": "2021-11-05T10:28:20Z",
		"timeModified": "2023-05-31T16:38:26Z"
	},
	"endDate": "2024-04-26T12:48:02Z",
	"id": "387b5178-732f-4706-9c41-e67a0cd00dc6",
	"ownerId": "f030a442-4aa3-3bf1-acf0-25f76194603f",
	"pipeline": {
		"analysisStorage": {
			"description": "1.2TB",
			"id": "6e1b6c8f-f913-48b2-9bd0-7fc13eda0fd0",
			"name": "Small",
			"ownerId": "8ec463f6-1acb-341b-b321-043c39d8716a",
			"tenantId": "f91bb1a0-c55f-4bce-8014-b2e60c0ec7d3",
			"tenantName": "ica-cp-admin",
			"timeCreated": "2021-11-05T10:28:20Z",
			"timeModified": "2023-05-31T16:38:26Z"
		},
		"code": "basic_pipeline",
		"description": "Reverses a fasta file and outputs to stdout.",
		"id": "bfecca03-6443-45bd-b313-e4f555cd0748",
		"language": "NEXTFLOW",
		"languageVersion": {
			"id": "2483549a-1530-4973-bb00-f3f6ccb7e610",
			"language": "NEXTFLOW",
			"name": "20.10.0"
		},
		"ownerId": "f030a442-4aa3-3bf1-acf0-25f76194603f",
		"pipelineTags": {
			"technicalTags": []
		},
		"tenantId": "02ea1bcf-6b20-4cbf-a9b2-724d1833eb07",
		"tenantName": "sbimb-wits",
		"timeCreated": "2024-04-26T12:23:56Z",
		"timeModified": "2024-04-26T12:23:56Z",
		"urn": "urn:ilmn:ica:pipeline:bfecca03-6443-45bd-b313-e4f555cd0748#basic_pipeline"
	},
	"reference": "regan_test_analysis_01-basic_pipeline-21ac67ed-ada1-4a82-8415-d2f83ec1e918",
	"startDate": "2024-04-26T12:34:19Z",
	"status": "SUCCEEDED",
	"summary": "",
	"tags": {
		"referenceTags": [],
		"technicalTags": [],
		"userTags": [
			"regan"
		]
	},
	"tenantId": "02ea1bcf-6b20-4cbf-a9b2-724d1833eb07",
	"tenantName": "sbimb-wits",
	"timeCreated": "2024-04-26T12:34:11Z",
	"timeModified": "2024-04-26T12:48:04Z",
	"userReference": "regan_test_analysis_01"
}
```
We can sort the project analysis by any one of the following fields:   
- reference
- userReference
- pipeline
- status
- startDate
- endDate
- summary   

We'll use `userReference` because this value is set by the user themselves.    

We'd like to use the CLI to trigger a pipeline run. In order to use the CLI to automate the running of pipelines, an _analysis code_ is required. To get this code, we need to run a pipeline to completion with the UI. We can then run the pipeline on any appropriate file (a small test file will be good enough). This is explained in more detail over [here](https://help.ica.illumina.com/tutorials/launchpipecli).    

Once this analysis code is extracted from the completed pipeline run, we can then use it in tandem with the _file id_ to execute the pipeline on a particular file. The command we use for a basic pipeline run (or analysis) is the following:
```bash
icav2 projectpipelines start nextflow $pipeline_id \
	--user-reference $user_reference \
	--project-id $project_id \
	--storage-size $storage_size \
	--input $file_ref
```
The `$file_ref` variable is constructed from the analysis code and the file id as `analysisCode:fileId`. The script for the more detailed process is [start_nextflow_pipeline_analysis.sh](bash/helper_scripts/start_nextflow_pipeline_analysis.sh).   