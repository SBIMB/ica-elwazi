## Delete Output File/Folder   
To delete a file from a project, the following CLI command can be used:
```bash
icav2 projectdata delete <path> --project-id <project_id>
```
The script [delete_file_by_path.sh](./../bash/helper_scripts/delete_file_by_path.sh) extracts the path of the file in the ICA storage and then proceeds to delete the file. Below is a screenshot of the process completed successfully:   
![Delete File from ICA Storage](./../public/assets/images/delete_file_by_path_script.png "Delete File from ICA Storage")  

Deleting a folder follows the same logic as deleting a file. For instance, the following command would delete a folder by its path in the ICA storage:
```bash
icav2 projectdata delete $folder_ica_storage_path --project-id $project_id
```
The script [delete_folder_by_path](bash/helper_scripts/delete_folder_by_path.sh) contains logic for extracting the folder path by filtering a list of project data by the folder name. The folder path is then used to delete the folder.    
