nextflow.enable.dsl=2

params.genes_bed        = null
params.biallelic_genotype_shards = null
params.anno_shards      = null
params.siteqc_shards    = null


workflow {
    
    row_indices = Channel.of(1,2,3,4,5,6)

    biallelic_ch, anno_ch, siteqc_ch, bed_ch = FIND_SHARDS(
        row_indices,
        file(params.genes_bed),
        file(params.biallelic_genotype_shards),
        file(params.anno_shards),
        file(params.siteqc_shards)
    )

shard_results = biallelic_ch
    .combine(anno_ch)
    .combine(siteqc_ch)
    .combine(bed_ch)

     staged_inputs = shard_results.map {  biallelic_genotype_shards, anno_shards, siteqc_shards, gene_bed_file ->

        tuple(
            biallelic_genotype_shards.readLines().collect { file(it.trim()) },
            anno_shards.readLines().collect { file(it.trim()) },
            siteqc_shards.readLines().collect { file(it.trim()) },
            gene_bed_file
        )
    }

    RUN_GENE(staged_inputs)
}

process FIND_SHARDS {
    input:
            val row_index
            path(gene_bed)
            path(biallelic_genotype_shards)
            path(anno_shards)
            path(siteqc_shards)
            

    output:
        path "*_biallelic_genotype_shards.txt"
        path "*_anno_shards.txt"
        path "*_siteqc_shards.txt"
        path "temp_*.bed"

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
        tuple
            path(biallelic_vcfs)
            path(anno_vcfs)
            path(siteqc_vcfs)
            path(genes_bed_file)

    output:
        path "*.tsv"

    script:
    """
        bash get_gene_specific_variants_AggV3.sh \
            ${biallelic_vcfs} \
            ${aanno_vcfs} \
            ${siteqc_vcfs}
            ${ggenes_bed_file}
    """
}
