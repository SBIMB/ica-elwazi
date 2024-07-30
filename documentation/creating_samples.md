## Creating Samples   
A sample allows us to group our resources together so that further actions can be performed in a more systematic manner. We can create a sample with the CLI using the following command:
```bash
icav2 projectsamples create $sample_name \
    --project-id $project_id \
    --description $sample_description \
    --user-tag $sample_user_tag \
    --technical-tag $sample_technical_tag
```
The script [create_project_sample.sh](./../bash/helper_scripts/create_project_sample.sh) contains the logic for creating a project sample. Once a sample is created, we can link it to an existing project (in addition to the project to which it already belongs). We can also upload files or folders to a sample.