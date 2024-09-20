#!/usr/bin/env nextflow

// using scanpy, perform automated scRNA analysis
process scanpy_analysis {
    container "dzezy/scanpy_benchmark"

    cache 'lenient'

    publishDir "$projectDir/results/${sample_name}/scanpy/", mode: 'copy'

    input:
    tuple val(sample_name), path(h5)

    output:
    path("scanpy_${sample_name}_benchmarks.log")
    path("scanpy_${sample_name}_plots.pdf")

    """
    time python ${projectDir}/scripts/scanpy_analysis.py ${h5} ${sample_name}
    """
}

// using seurat, perform automated scRNA analysis
process seurat_analysis {
    container "dzezy/seurat_benchmark"

    cache 'lenient'

    publishDir "$projectDir/results/${sample_name}/seurat/", mode: 'copy'

    input:
    tuple val(sample_name), path(h5)

    output:
    path("seurat_${sample_name}_benchmarks.log")
    path("seurat_${sample_name}_plots.pdf")

    """
    time Rscript ${projectDir}/scripts/seurat_analysis.R ${h5} ${sample_name}
    """
}

// define the workflow steps
workflow {
    h5_files = Channel.fromPath(params.h5_files)
    h5_files
        .map{ it -> [it.getBaseName(), it] }
        .set {named_h5}
    
    scanpy_analysis(named_h5)
    seurat_analysis(named_h5)

}
