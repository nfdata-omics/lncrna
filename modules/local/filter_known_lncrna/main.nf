process FILTER_KNOWN_LNCRNA {
tag "$meta.id"
label 'process_medium'

container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/gffcompare:0.12.6--h9948957_4' :
    'biocontainers/gffcompare:0.12.10--h9948957_0' }"

input:
tuple val(meta), path(novel_lncrna_gtf)
path known_lncrna_gtf

output:
tuple val(meta), path("*.truly_novel.gtf"), emit: gtf
path "*.tmap"                              , emit: tmap
path "*.stats"                             , emit: stats
path "versions.yml"                        , emit: versions

when:
task.ext.when == null || task.ext.when

script:
def prefix = task.ext.prefix ?: "${meta.id}"
"""
# Compare with known lncRNAs using gffcompare
gffcompare \\
    -G \\
    -o filter \\
    -r $known_lncrna_gtf \\
    -p $task.cpus \\
    $novel_lncrna_gtf

# Extract only truly novel lncRNAs (class codes 'u' or 'x')
# u = unknown intergenic transcript
# x = exonic overlap with reference on opposite strand
awk '\$3 == "u" || \$3 == "x" {print \$5}' filter.${novel_lncrna_gtf}.tmap | \\
    sort | uniq | \\
    grep -Ff - $novel_lncrna_gtf > ${prefix}.truly_novel.gtf

# Copy outputs for tracking
cp filter.${novel_lncrna_gtf}.tmap ${prefix}.tmap
cp filter.stats ${prefix}.stats

cat <<-END_VERSIONS > versions.yml
"${task.process}":
    gffcompare: \$(gffcompare --version 2>&1 | sed 's/gffcompare v//')
END_VERSIONS
"""
}
