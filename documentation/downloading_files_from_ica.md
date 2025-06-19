## Downloading Files from Illumina Connected Analytics   
The output data files can be downloaded from the ICA storage using the web UI. A download can even be scheduled through the web UI.   

![Schedule Download through Web UI](./../public/assets/images/download_scheduled_with_ui.png "Schedule Download through Web UI")   

As mentioned above, we would like to trigger a download process as soon as the status reaches **"SUCCEEDED"**. After a successful download, we can then delete the files from the ICA storage.   

Files can be downloaded from the ICA storage using the CLI with the following CLI command:
```bash
icav2 projectdata list # take note of the sourcePath
icav2 projectdata download <sourcePath or dataId> <targetPath>
```
The script [download_file_by_path.sh](./../bash_scripts/download_file_by_path.sh) can be tested with some test data. A successful implementation of this download script should look as follows:   

![Download File from ICA Storage](./../public/assets/images/successful_download_script.png "Download File from ICA Storage")  