#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

process wait_for_ica {
    input:
    path config
    each path(analysis_file)

    output:
    env 'status'

    script:
    """
      project_id=\$(cat ${config} | jq -r .ica_job_project.id)
      job_id=\$(cat ${analysis_file} | jq -r .id)

      timeout=\$((24*60*60)) # one day
      start=\$SECONDS
      request_errors=0

      while true; do
        elapsed=\$(( SECONDS - start ))

        if [ \$elapsed -ge \$timeout ]; then
          echo "Job \${job_id} timed out after \$timeout seconds."
          exit 1
        fi

        if [ \$request_errors -ge 10 ]; then
          echo "Failing \${job_id} after 10 failed status requests."
          exit 1
        fi        

        # Assume a transient network error
        if ! icav2 projectanalyses get --output-format json --project-id \${project_id} \${job_id} 1> job.json ; then
            request_errors=\$(( request_errors + 1 ))
            sleep \$(( 60 * request_errors ))
            continue
        fi
        request_errors=0

        status=\$(cat job.json | jq -r .status )

        if [[ \${status} == SUCCEEDED ]]; then
          break
        fi

        if [[ \${status} == FAILED ]]; then
          echo "Job \${job_id} failed."
          cat job.json | jq -r .summary
          exit 1
        fi

        sleep 60
      done
    """
}
