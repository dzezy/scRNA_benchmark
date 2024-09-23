#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// using scanpy, perform automated scRNA analysis
process scanpy_analysis {
    container "dzezy/scanpy_benchmark"

    cache 'lenient'

    publishDir "${params.outputDir}/${sample_name}/benchmarks/", mode: 'copy', pattern: '*.log'
    publishDir "${params.outputDir}/${sample_name}/plots/", mode: 'copy', pattern: '*.pdf'

    input:
    tuple val(sample_name), path(h5)

    output:
    path("scanpy_${sample_name}_benchmarks.log")
    path("scanpy_${sample_name}_plots.pdf")

    """
    time python ${params.scriptDir}/scanpy_analysis.py ${h5} ${sample_name}
    """
}

// using seurat, perform automated scRNA analysis
process seurat_analysis {
    container "dzezy/seurat_benchmark"
    
    cache 'lenient'

    publishDir "${params.outputDir}/${sample_name}/benchmarks/", mode: 'copy', pattern: '*.log'
    publishDir "${params.outputDir}/${sample_name}/plots/", mode: 'copy', pattern: '*.pdf'

    input:
    tuple val(sample_name), path(h5)

    output:
    path("seurat_${sample_name}_benchmarks.log")
    path("seurat_${sample_name}_plots.pdf")

    """
    time Rscript ${params.scriptDir}/seurat_analysis.R ${h5} ${sample_name}
    """
}

// define the workflow steps
workflow {
    // extract h5_files and their basename
    h5_files = Channel.fromPath(params.h5_files)
    h5_files
        .map{ it -> [it.getBaseName(), it] }
        .set {named_h5}
    
    scanpy_analysis(named_h5)
    seurat_analysis(named_h5)

}
