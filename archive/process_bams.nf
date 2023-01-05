nextflow.enable.dsl=2
params.bam_dir = '/mnt/c/desktop/ivar_bam'
// scripts are in scripts folder relative to the process_bams.nf file
script_dir = file("$baseDir/scripts")

// get bam and bai files
bamList = Channel.fromPath("${params.bam_dir}/*.bam")
bamIndex = Channel.fromPath("${params.bam_dir}/*.bam.bai")





process bamToMutations {
    publishDir "results", mode: 'copy'

    input:
    file bam
    file bamIndex

    output:
    file "${bam}.tsv" 

// scripts/extract_bases
    script:
    """
    python ${script_dir}/extract_bases.py -i ${bam} -o ${bam}.tsv
    """
}

workflow {
    bamToMutations(bamList, bamIndex)
   
}