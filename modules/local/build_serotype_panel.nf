process BUILD_SEROTYPE_PANEL {
    tag "${reference_panel.name}"
    label 'process_low'
    conda "bioconda::minimap2=2.26"
    container "quay.io/biocontainers/minimap2:2.26--he4a0461_2"

    input:
    path reference_dir
    path reference_panel

    output:
    path "denv_ref.mmi",                  emit: index
    path "all_denv_references_tagged.fa", emit: tagged_fasta
    path "versions.yml",                  emit: versions

    script:
    """
    : > all_denv_references_tagged.fa
    tail -n +2 ${reference_panel} | while IFS=\$'\\t' read -r acc fasta sero sylvatic ds desc; do
        [ -z "\$acc" ] && continue
        fpath="${reference_dir}/\$fasta"
        if [ ! -s "\$fpath" ]; then
            echo "ERROR: reference file missing for \$acc: \$fpath" >&2
            exit 1
        fi
        tag="\${sero}"
        if [ "\$sylvatic" = "true" ]; then
            tag="\${sero}_SYLVATIC"
        fi
        awk -v hdr=">\${tag}|\${acc}" 'NR==1 {print hdr; next} {print}' "\$fpath" >> all_denv_references_tagged.fa
    done

    if [ ! -s all_denv_references_tagged.fa ]; then
        echo "ERROR: combined reference panel is empty" >&2
        exit 1
    fi

    minimap2 -d denv_ref.mmi all_denv_references_tagged.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version)
    END_VERSIONS
    """

    stub:
    """
    touch denv_ref.mmi all_denv_references_tagged.fa
    printf '"%s":\\n    minimap2: stub\\n' "${task.process}" > versions.yml
    """
}
