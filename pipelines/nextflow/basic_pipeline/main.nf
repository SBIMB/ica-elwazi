#!/usr/bin/env nextflow

params.in = "$HOME/sample.fa"

sequences = file(params.in)
SPLIT = (System.properties['os.name'] == 'macOS' ? 'gcsplit' : 'csplit')

process splitSequences {

    container 'public.ecr.aws/lts/ubuntu:22.04'

    input:
    file 'input.fa' from sequences

    output:
    file 'seq_*' into records

    """
    $SPLIT input.fa '%^>%' '/^>/' '{*}' -f seq_
    """

}

process reverse {
    
    container 'public.ecr.aws/lts/ubuntu:22.04'
    publishDir 'out'

    input:
    file x from records
    
    output:
    file 'test.txt'

    """
    cat $x | rev > test.txt
    """
}