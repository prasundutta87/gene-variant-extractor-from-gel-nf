nextflow.enable.dsl=2

params.genes_bed                  = null
params.biallelic_genotype_shards  = null
params.anno_shards                = null
params.siteqc_shards              = null


process FIND_SHARDS {

    input:
        val row_index
        path genes_bed
        path biallelic_genotype_shards
        path anno_shards
        path siteqc_shards

    output:
        tuple path("*_biallelic_genotype_shards.txt"),
              path("*_anno_shards.txt"),
              path("*_siteqc_shards.txt"),
              path("temp_*.bed")

    script:
    """
    bash find_gene_shards.sh \
        ${row_index} \
        ${genes_bed} \
        ${biallelic_genotype_shards} \
        ${anno_shards} \
        ${siteqc_shards}
    """
}


process RUN_GENE {

    container "prasundutta87/gene-variant-extractor-from-gel-docker-image:1.0.0"

    publishDir "results", mode: 'copy'

    input:
        tuple path(biallelic_txt),
              path(anno_txt),
              path(siteqc_txt),
              path(gene_bed_file)

    output:
        path "*.tsv"

    script:
    """
    bash get_gene_specific_variants_AggV3.sh \
        ${biallelic_txt} \
        ${anno_txt} \
        ${siteqc_txt} \
        ${gene_bed_file}
    """
}


workflow {

    row_indices = Channel.of(6)

    shard_results = FIND_SHARDS(
        row_indices,
        file(params.genes_bed),
        file(params.biallelic_genotype_shards),
        file(params.anno_shards),
        file(params.siteqc_shards)
    )

    RUN_GENE(shard_results)
}