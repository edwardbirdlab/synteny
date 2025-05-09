process JCVI {

    label 'process_single'
    tag "$sample_id"
    container = 'quay.io/ecoflowucl/jcvi:python-3.10_last-1522'

    input:
    tuple val(sample_id), path(fasta), path(gff)

    output:
    tuple val(sample_id), path( "${sample_id}.cds" ), path( "${sample_id}.bed" ) , emit: new_format
    path( "${sample_id}.bed" ) , emit: beds
    path "versions.yml", emit: versions

    script:
    """
    #Run the basic transformation of gff to bed and fasta to cds conversions.
    python -m jcvi.formats.gff bed --type=${params.jcvi_bed_type} --key=ID ${gff} -o ${sample_id}.bed
    python -m jcvi.formats.fasta format ${fasta} ${sample_id}.cds

    md5sum "${sample_id}.bed" > "${sample_id}.bed.md5"
    md5sum "${sample_id}.cds" > "${sample_id}.cds.md5"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Python version: \$(python --version  | sed 's/[^0-9]*//')
        JCVI \$(pip show jcvi | grep "Version:")
    END_VERSIONS
    """
}
