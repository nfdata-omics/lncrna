process MERGE_FINAL_ANNOTATION {
tag "$meta.id"
label 'process_low'

container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gffread:0.12.7--h077b44d_5' :
        'biocontainers/gffread:0.12.7--h077b44d_6' }"

input:
tuple val(meta), path(lncrna_gtf)       // RENAME_LNCRNA.out.gtf
path protein_coding_gtf                 // GTF_FILTER_PROTEIN_CODING.out.protein_gtf
path fasta                              // PREPARE_GENOME.out.fasta

output:
tuple val(meta), path("*.final_all.gtf")        , emit: gtf
tuple val(meta), path("*.final_all.fa")         , emit: fasta
tuple val(meta), path("*.lncrna_only.fa")       , emit: lncrna_fasta
tuple val(meta), path("*.protein_only.fa")      , emit: protein_fasta
path "versions.yml"                             , emit: versions

when:
task.ext.when == null || task.ext.when

script:
def prefix = task.ext.prefix ?: "${meta.id}"
"""
# Merge lncRNA and protein-coding GTFs
cat $lncrna_gtf $protein_coding_gtf > ${prefix}.final_all.gtf

# Extract sequences for merged annotation
gffread ${prefix}.final_all.gtf \\
    -g $fasta \\
    -w ${prefix}.final_all.fa \\
    -W

# Extract lncRNA sequences only
gffread $lncrna_gtf \\
    -g $fasta \\
    -w ${prefix}.lncrna_only.fa \\
    -W

# Extract protein-coding sequences only
gffread $protein_coding_gtf \\
    -g $fasta \\
    -w ${prefix}.protein_only.fa \\
    -W

cat <<-END_VERSIONS > versions.yml
"${task.process}":
    gffread: \$(gffread --version 2>&1)
END_VERSIONS
"""
}
