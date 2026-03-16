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
        path "*_biallelic_genotype_shards.txt", emit: biallelic
        path "*_anno_shards.txt", emit: anno
        path "*_siteqc_shards.txt", emit: siteqc
        path "temp_*.bed", emit: gene_bed

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
        path(biallelic_txt)
        path(anno_txt)
        path(siteqc_txt)
        path(gene_bed_file)
        path(staged_shards)

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

    FIND_SHARDS(
        row_indices,
        file(params.genes_bed),
        file(params.biallelic_genotype_shards),
        file(params.anno_shards),
        file(params.siteqc_shards)
    )

    ch_biallelic = FIND_SHARDS.out.biallelic
        .ifEmpty { exit 1, "Cannot find file : ${FIND_SHARDS.out.biallelic}" }
        .splitText { it.trim() }
        .filter { it != "" }
        .map { file(it) }

    ch_anno = FIND_SHARDS.out.anno
        .ifEmpty { exit 1, "Cannot find file : ${FIND_SHARDS.out.anno}" }
        .splitText { it.trim() }
        .filter { it != "" }
        .map { file(it) }

    ch_siteqc = FIND_SHARDS.out.siteqc
        .ifEmpty { exit 1, "Cannot find file : ${FIND_SHARDS.out.siteqc}" }
        .splitText { it.trim() }
        .filter { it != "" }
        .map { file(it) }

    shard_results = ch_biallelic.concat(ch_anno).concat(ch_siteqc)

    RUN_GENE(
        FIND_SHARDS.out.biallelic,
        FIND_SHARDS.out.anno,
        FIND_SHARDS.out.siteqc,
        FIND_SHARDS.out.gene_bed,
        shard_results.unique().collect()
    )
}