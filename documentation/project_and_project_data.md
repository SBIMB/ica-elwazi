## Project and Project Data   
A project can be created in the UI. After a project is created, the project object can be obtained by using the following CLI command:
```bash
projects=$(icav2 projects list)
```
The response body will contain a JSON object with a list of existing projects. A single project object contains a lot of details about the project. (The object is too big to be displayed here).   

To get the details of a specific project, the project id or name will be required. For example, we have a [script](bash/helper_scripts/get_project_id.sh) in the `/bash` directory that extracts the `project_id` by using the `project_name` (which we create when creating a project).   

We make use of the [jq](https://jqlang.github.io/jq/) library (which is a JSON parser) in several of our scripts.  