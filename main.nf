nextflow.enable.dsl=2

params.genes_bed                  = null
params.biallelic_genotype_shards  = null
params.anno_shards                = null
params.siteqc_shards              = null
params.sample_list                = null
params.seq_report_100k            = null
params.seq_report_gms             = null
params.rd_phenotype               = null
params.gms_phenotype              = null
params.outdir                     = "results"


process FIND_SHARDS {
    input:
        val  row_index
        path genes_bed
        path biallelic_genotype_shards
        path anno_shards
        path siteqc_shards

    output:
        path "*_biallelic_genotype_shards.txt",     emit: biallelic
        path "*_anno_shards.txt",                   emit: anno
        path "*_siteqc_shards.txt",                 emit: siteqc
        path "*_biallelic_genotype_shards_idx.txt", emit: biallelic_idx
        path "*_anno_shards_idx.txt",               emit: anno_idx
        path "*_siteqc_shards_idx.txt",             emit: siteqc_idx
        path "temp_*.bed",                          emit: gene_bed

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
    container "prasundutta87/gene-variant-extractor-from-gel-docker-image:2.0.0"

    stageInMode 'copy'

    publishDir "${params.outdir}", mode: 'copy', pattern: "*.tsv"

    input:
        path gene_bed_file
        // stageAs gives each file a unique name using its index
        // VCF and its index are staged with matching basenames so bcftools finds them
        path biallelic_vcfs, stageAs: "biallelic_??.vcf.gz"
        path biallelic_idx,  stageAs: "biallelic_??.vcf.gz.tbi"
        path anno_vcfs,      stageAs: "anno_??.vcf.gz"
        path anno_idx,       stageAs: "anno_??.vcf.gz.tbi"
        path siteqc_vcfs,    stageAs: "siteqc_??.vcf.gz"
        path siteqc_idx,     stageAs: "siteqc_??.vcf.gz.tbi"

    output:
        tuple path(gene_bed_file), path("*.tsv"), emit: gene_tsv

    script:
    """
    # Write staged filenames to text files for the bash script
    ls biallelic_??.vcf.gz > biallelic_staged.txt
    ls anno_??.vcf.gz      > annotation_staged.txt
    ls siteqc_??.vcf.gz    > siteqc_staged.txt

    bash get_gene_specific_variants_AggV3.sh \
        biallelic_staged.txt \
        annotation_staged.txt \
        siteqc_staged.txt \
        ${gene_bed_file}
    """
}


process ANNOTATE_GENE {
    container "prasundutta87/gene-variant-extractor-from-gel-docker-image:2.0.0"

    publishDir "${params.outdir}", mode: 'copy'

    input:
        path gene_tsv
        path sample_list
        path seq_report_100k
        path seq_report_gms
        path rd_phenotype
        path gms_phenotype

    output:
        path "*/*_for_review.tsv"

    script:
    def gene_name = gene_tsv.baseName
    """
    join_additional_anno_with_genes \
        ${gene_name} \
        ${sample_list} \
        ${seq_report_100k} \
        ${seq_report_gms} \
        ${rd_phenotype} \
        ${gms_phenotype}
    """
}


// Reads shard path text file and returns [names, files] tuple for staging
def processShardFiles(channel) {
    return channel
        .splitText { it.trim() }
        .filter   { it != "" }
        .map      { line -> file(line.trim()) }
        .collect()
}


workflow {

    def n_genes = file(params.genes_bed).countLines()
    ch_row_indices = Channel.from( 1..n_genes )

    FIND_SHARDS(
        ch_row_indices,
        file(params.genes_bed),
        file(params.biallelic_genotype_shards),
        file(params.anno_shards),
        file(params.siteqc_shards)
    )

    // Stage VCF and index files via Nextflow using Lifebit S3 credentials
    ch_biallelic     = processShardFiles(FIND_SHARDS.out.biallelic)
    ch_anno          = processShardFiles(FIND_SHARDS.out.anno)
    ch_siteqc        = processShardFiles(FIND_SHARDS.out.siteqc)
    ch_biallelic_idx = processShardFiles(FIND_SHARDS.out.biallelic_idx)
    ch_anno_idx      = processShardFiles(FIND_SHARDS.out.anno_idx)
    ch_siteqc_idx    = processShardFiles(FIND_SHARDS.out.siteqc_idx)

    RUN_GENE(
        FIND_SHARDS.out.gene_bed,
        ch_biallelic,
        ch_anno,
        ch_siteqc,
        ch_biallelic_idx,
        ch_anno_idx,
        ch_siteqc_idx
    )

    ANNOTATE_GENE(
        RUN_GENE.out.gene_tsv.map { bed, tsv -> tsv },
        file(params.sample_list),
        file(params.seq_report_100k),
        file(params.seq_report_gms),
        file(params.rd_phenotype),
        file(params.gms_phenotype)
    )
}


workflow.onComplete {
    log.info "Pipeline completed : ${workflow.complete}"
    log.info "Duration           : ${workflow.duration}"
    log.info "Success            : ${workflow.success}"
    log.info "Output directory   : ${params.outdir}"
}
