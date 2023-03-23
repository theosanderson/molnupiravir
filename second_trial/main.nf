#!/usr/bin/env nextflow

params.reads = "seqs/*_{1,2}.fastq.gz"
params.ref = "ref.fa"
params.outdir = "results"
params.outdir_tsv = "results_tsv"

// Channel for input reads
Channel
    .fromFilePairs(params.reads, size: 2)
    .set { read_pairs }

// Channel for the reference genome
Channel
    .fromPath(params.ref)
    .set { ref_genome }

// Process to index the reference genome
process index_ref {
    input:
    path ref from ref_genome

    output:
    path "${ref}.mmi" into indexed_ref

    script:
    """
    minimap2 -d "${ref}.mmi" "${ref}"
    """
}

// Process to map reads using minimap2
process map_reads {
    input:
    tuple val(sample), path(reads) from read_pairs
    path indexed_ref from indexed_ref.collect()

    output:
    path "${sample}.sam" into mapped_reads

    script:
    """
    minimap2 -ax sr "${indexed_ref}" "${reads[0]}" "${reads[1]}" > "${sample}.sam"
    """
}

// Process to convert SAM to sorted BAM and store results
process store_results {
    publishDir params.outdir, mode: 'copy'

    input:
    path mapped_sam from mapped_reads



    output:
    path "${mapped_sam.simpleName}.sorted.ba*" into final_results; 

    script:
    """
    samtools view -S -b "${mapped_sam}" | samtools sort -o "${mapped_sam.simpleName}.sorted.bam"
    samtools index "${mapped_sam.simpleName}.sorted.bam"
    """
}


// Process to extract bases using extract_bases.py. We need to just pass the BAM file to the script, not the BAI
process extract_bases {
    publishDir params.outdir_tsv, mode: 'copy'

    input:
    path bam_file from final_results
    
   
    output:
    path "*.tsv" into extracted_bases

    script:
    """
    # we need to jsut get the bam

    python ${baseDir}/extract_bases.py -i "${bam_file[0]}" -o "${bam_file[0]}.tsv"
    """
}