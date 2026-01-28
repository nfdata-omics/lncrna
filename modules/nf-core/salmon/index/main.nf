process SALMON_INDEX {
    tag "$transcript_fasta"
    label "process_medium"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/salmon:1.10.3--h6dccd9a_2' :
        'biocontainers/salmon:1.10.3--h6dccd9a_2' }"

    input:
    path genome_fasta
    path transcript_fasta

    output:
    path "index"       , emit: index
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def decoys = ''
    def fasta = transcript_fasta
    if (genome_fasta) {
        if ("${genome_fasta}".endsWith('.gz')) {
            genome_fasta = "<(gunzip -c ${genome_fasta})"
        }
        decoys='-d decoys.txt'
        fasta='gentrome.fa'
    }
    if ("${transcript_fasta}".endsWith('.gz')) {
        transcript_fasta = "<(gunzip -c ${transcript_fasta})"
    }
    """
    if [ -n '$genome_fasta' ]; then
        grep '^>' $genome_fasta | cut -d ' ' -f 1 | cut -d \$'\\t' -f 1 | sed 's/>//g' > decoys.txt
        cat $transcript_fasta $genome_fasta > $fasta
    fi

    salmon \\
        index \\
        --threads $task.cpus \\
        -t $fasta \\
        $decoys \\
        $args \\
        -i index

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        salmon: \$(echo \$(salmon --version) | sed -e "s/salmon //g")
    END_VERSIONS
    """

    stub:
    """
    mkdir index
    touch index/complete_ref_lens.bin
    touch index/ctable.bin
    touch index/ctg_offsets.bin
    touch index/duplicate_clusters.tsv
    touch index/info.json
    touch index/mphf.bin
    touch index/pos.bin
    touch index/pre_indexing.log
    touch index/rank.bin
    touch index/refAccumLengths.bin
    touch index/ref_indexing.log
    touch index/reflengths.bin
    touch index/refseq.bin
    touch index/seq.bin
    touch index/versionInfo.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        salmon: \$(echo \$(salmon --version) | sed -e "s/salmon //g")
    END_VERSIONS
    """
}
