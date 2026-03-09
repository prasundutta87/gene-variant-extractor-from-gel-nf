nextflow.enable.dsl=2

params.genes_bed        = null
params.biallelic_shards = null
params.anno_shards      = null
params.siteqc_shards    = null


workflow {
    
    row_indices = Channel.of(1,2,3,4,5,6)

    RUN_GENE(
        row_indices,
        file(params.genes_bed),
        file(params.biallelic_shards),
        file(params.anno_shards),
        file(params.siteqc_shards)
    )
}

process RUN_GENE {

    container "prasundutta87/gene-variant-extractor-from-gel-docker-image:1.0.0"

    publishDir "results", mode: 'copy'

    input:
        val row_index
        path genes_bed
        path biallelic_shards
        path anno_shards
        path siteqc_shards

    output:
        path "*.tsv"

    script:
    """
    bash ${projectDir}/scripts/get_gene_specific_variants_AggV3.sh \
        ${row_index} \
        ${genes_bed} \
        ${biallelic_shards} \
        ${anno_shards} \
        ${siteqc_shards}
    """
}
