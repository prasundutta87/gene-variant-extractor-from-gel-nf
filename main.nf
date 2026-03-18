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
        path "*_biallelic_genotype_shards_idx.txt", emit: biallelic_idx
        path "*_anno_shards_idx.txt", emit: anno_idx
        path "*_siteqc_shards_idx.txt", emit: siteqc_idx

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
        path(gene_bed_file)
        tuple val(biallelic_names), path('biallelic_*')
        tuple val(anno_names), path('annotation_*')
        tuple val(siteqc_names), path('siteqc_*')
        tuple val(biallelic_idx_names), path('biallelic_idx_*')
        tuple val(anno_idx_names), path('annotation_idx_*')
        tuple val(siteqc_idx_names), path('siteqc_idx_*')
        path combine_duckplyr

    output:
        path "*.tsv"

    script:
    // Rename staged files to their proper names  
    def rename_cmds_b = biallelic_names.withIndex().collect { name, idx ->
        "mv biallelic_${idx+1} ${name}"
    }.join('\n    ')
    def rename_cmds_a = anno_names.withIndex().collect { name, idx ->
        "mv annotation_${idx+1} ${name}"
    }.join('\n    ')
    def rename_cmds_s = siteqc_names.withIndex().collect { name, idx ->
        "mv siteqc_${idx+1} ${name}"
    }.join('\n    ')
    // rename indexes
    def rename_cmds_b_idx = biallelic_idx_names.withIndex().collect { name, idx ->
        "mv biallelic_idx_${idx+1} ${name}"
    }.join('\n    ')
    def rename_cmds_a_idx = anno_idx_names.withIndex().collect { name, idx ->
        "mv annotation_idx_${idx+1} ${name}"
    }.join('\n    ')
    def rename_cmds_s_idx = siteqc_idx_names.withIndex().collect { name, idx ->
        "mv siteqc_idx_${idx+1} ${name}"
    }.join('\n    ')
    // Create new text files with staged filenames
    def biallelic_list = biallelic_names.collect { "${it}" }.join('\n')
    def anno_list = anno_names.collect { "${it}" }.join('\n')
    def siteqc_list = siteqc_names.collect { "${it}" }.join('\n')
    """
    # Rename files to their unique names
    ${rename_cmds_b}
    ${rename_cmds_a}
    ${rename_cmds_s}
    ${rename_cmds_b_idx}
    ${rename_cmds_a_idx}
    ${rename_cmds_s_idx}

    # Create new text files pointing to staged filenames
    echo "${biallelic_list}" > biallelic_staged.txt
    echo "${anno_list}" > annotation_staged.txt
    echo "${siteqc_list}" > siteqc_staged.txt

    bash get_gene_specific_variants_AggV3.sh \
        biallelic_staged.txt \
        annotation_staged.txt \
        siteqc_staged.txt \
        ${gene_bed_file}
    """
}

// Reusable function to process shard files into unique named tuples
def processShardFiles(channel, channelName) {
    return channel
        .ifEmpty { exit 1, "Cannot find file : ${channelName}" }
        .splitText { it.trim() }
        .filter { it != "" }
        .map { line ->
            // extract s3://... or file path token
            def m = (line =~ /(s3:\/\/\S+)/)
            def uri = m ? m[0][1] : line.trim()
            // remove the s3:// and split to get parent dir and base name
            def pathNoSchema = uri.replaceFirst(/^s3:\/\//,'')
            def parts = pathNoSchema.tokenize('/')
            // join all parts with "_"
            def newName = parts.join('_')
            // emit tuple (name, uri)
            return [ newName, uri ]
        }
        // dedupe by the newName so we only keep one row per generated name
        .unique { it[0] }
        // Group into a single emission with list of names and list of file objects
        .toList()
        .map { tuples -> 
            def names = tuples.collect { it[0] }
            def files = tuples.collect { file(it[1]) }
            println "Collected ${names.size()} unique ${channelName} files"
            [ names, files ]
        }
}


workflow {

    row_indices = Channel.of(6)
    combine_duckplyr = Channel.fromPath("${workflow.projectDir}/bin/combine_variant_using_duckplyr.R")

    FIND_SHARDS(
        row_indices,
        file(params.genes_bed),
        file(params.biallelic_genotype_shards),
        file(params.anno_shards),
        file(params.siteqc_shards)
    )

    ch_biallelic = processShardFiles(FIND_SHARDS.out.biallelic, 'biallelic')
    ch_annot = processShardFiles(FIND_SHARDS.out.anno, 'annotation')
    ch_siteqc = processShardFiles(FIND_SHARDS.out.siteqc, 'siteqc')
    ch_biallelic_idx = processShardFiles(FIND_SHARDS.out.biallelic_idx, 'biallelic_idx')
    ch_annot_idx = processShardFiles(FIND_SHARDS.out.anno_idx, 'annotation_idx')
    ch_siteqc_idx = processShardFiles(FIND_SHARDS.out.siteqc_idx, 'siteqc_idx')

    RUN_GENE(
        FIND_SHARDS.out.gene_bed,
        ch_biallelic,
        ch_annot,
        ch_siteqc,
        ch_biallelic_idx,
        ch_annot_idx,
        ch_siteqc_idx,
        combine_duckplyr
    )
}