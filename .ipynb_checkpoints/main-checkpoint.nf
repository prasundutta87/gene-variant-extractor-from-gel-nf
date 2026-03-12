nextflow.enable.dsl=2

params.genes_bed        = null
params.biallelic_genotype_shards = null
params.anno_shards      = null
params.siteqc_shards    = null

process FIND_SHARDS {
    input:
            val row_index
            path(gene_bed)
            path(biallelic_genotype_shards)
            path(anno_shards)
            path(siteqc_shards)
            

   output:
    path "*_biallelic_genotype_shards.txt", emit: biallelic
    path "*_anno_shards.txt", emit: anno
    path "*_siteqc_shards.txt", emit: siteqc
    path "temp_*.bed", emit: bed

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

workflow {
    
    row_indices = Channel.of(1,2,3,4,5,6)

    shard_proc = FIND_SHARDS(
    row_indices,
    file(params.genes_bed),
    file(params.biallelic_genotype_shards),
    file(params.anno_shards),
    file(params.siteqc_shards)
)

biallelic_ch = shard_proc.out.biallelic
anno_ch      = shard_proc.out.anno
siteqc_ch    = shard_proc.out.siteqc
bed_ch       = shard_proc.out.bed

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