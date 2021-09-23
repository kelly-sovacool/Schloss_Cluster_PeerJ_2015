
rule download_vsearch_altVersion:
    output:
        tar='code/vsearch-{version}-linux-x86_64.tar.gz',
        bin='code/vsearch-{version}/vsearch'
    params:
        bin='bin/vsearch-{version}-linux-x86_64/bin/vsearch'
    shell:
        """
        wget -P bin/ https://github.com/torognes/vsearch/releases/download/v{wildcards.version}/vsearch-{wildcards.version}-linux-x86_64.tar.gz
        tar -C bin/ -xzvf {output.tar}
        mv {params.bin} {output.bin}
        rm -rf bin/
        """

rule get_vdgc_sensspec:
    output:
        sensspec="data/miseq/miseq_0.2_01.vdgc.sensspec"
    shell:
        """
        make -j 8 {output.sensspec}
        """