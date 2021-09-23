
rule download_vsearch_altVersion:
    output:
        bin='code/vsearch-{version}/vsearch'
    params:
        tar='bin/vsearch-{version}-linux-x86_64.tar.gz',
        bin='bin/vsearch-{version}-linux-x86_64/bin/vsearch'
    shell:
        """
        wget -P bin/ https://github.com/torognes/vsearch/releases/download/v{wildcards.version}/vsearch-{wildcards.version}-linux-x86_64.tar.gz
        tar -C bin/ -xzvf {params.tar}
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

rule remove_gaps_query:
    input:
        fna="data/miseq_vdgc_debug/miseq_1.0_01/{dataset}.fasta" # TODO
    output:
        degap=temp("data/miseq_vdgc_debug/miseq_1.0_01.tmp.fasta"),
        fna="data/miseq_vdgc_debug/miseq_1.0_01.ng.fasta"
    log:
        'log/miseq_vdgc_debug/miseq_1.0_01/remove_gaps_query.log'
    params:
        outdir="data/miseq_vdgc_debug/miseq_1.0_01/",
        prefix='miseq_1.0_01.tmp'
    resources:
        procs=2
    shell:
        """
        mothur '#set.dir(output={params.outdir}); set.logfile(name={log});
            degap.seqs(fasta={input.fna}, processors={resources.procs});
            rename.file(fasta=current, prefix={params.prefix}) '
        # vsearch doesn't support dots or hyphens in fasta headers
        cat {output.degap} | sed 's/[\.-]/_/g' > {output.fna}
        """


rule vsearch_sort:
    input:
        fna="data/miseq_vdgc_debug/miseq_1.0_01.ng.fasta" #data/miseq/miseq_1.0_01.ng.fasta
    output:
        fna="data/miseq_vdgc_debug/miseq_1.0_01.ng.sorted.fasta",
        uc="data/miseq_vdgc_debug/miseq_1.0_01.ng.sorted.uc"
    shell:
        """
        vsearch \
            --derep_fulllength {input.fna} \
            --sizeout \
            --minseqlength 30 \
            --threads 1 \
            --uc {output.uc} \
            --output {output.fna} \
            --strand both
        """

rule vsearch_de_novo:
    input:
        query=rules.vsearch_sort.output.fna
    output:
        uc='results/miseq_1.0_01/de_novo/miseq_1.0_01.uc'
    benchmark:
        'benchmarks/miseq_1.0_01/vsearch.method_de_novo.miseq_1.0_01.txt'
    params:
        perc_identity=perc_identity,
        min_seq_length=min_seq_length,
        max_accepts=max_accepts,
        max_rejects=max_rejects,
        word_length=word_length
    resources:
        procs=8
    shell:
        """
        vsearch --cluster_smallmem {input.query} \
            --usersort \
            --uc {output.uc} \
            --threads {resources.procs} \
            --id {params.perc_identity} \
            --minseqlength {params.min_seq_length} \
            --maxaccepts {params.max_accepts} \
            --maxrejects {params.max_rejects} \
            --wordlength {params.word_length} \
            --strand both \
            --relabel OTU_
        """

rule uc_to_list:
    input:
        code='code/R/uc_to_list.R',
        sorted=rules.vsearch_sort.output.uc,
        clustered='results/miseq_1.0_01/{method}/{dataset_step}.uc'
    output:
        list='results/miseq_1.0_01/{method}/{dataset_step}.list'
    script:
        'code/R/uc_to_list.R'

# vsearch doesn't support dots or hyphens in sequence names.
rule prep_count_table:
    input:
        count_table=prep_samples("data/{dataset}/processed/{dataset}.count_table") # TODO
    output:
        count_table='data/miseq_vdgc_debug/miseq_1.0_01.count_table'
    shell:
        """
        cat {input.count_table} |  sed 's/[\.-]/_/g' > {output.count_table}
        """

# replace dots/hyphens in fasta headers (1st & 2nd occurrence),
# but not in the distance value (3rd occurrence).
rule prep_dist:
    input:
        dist=prep_samples("results/{dataset}/{dataset}.dist") # TODO
    output:
        dist='data/miseq_1.0_01/miseq_1.0_01.dist'
    shell:
        """
        cat {input.dist} |  sed 's/[\.-]/_/' | sed 's/[\.-]/_/' > {output.dist}
        """

rule sensspec_vsearch:
    input:
        list="results/miseq_1.0_01/{method}/miseq_1.0_01.list",
        count_table='data/miseq_vdgc_debug/miseq_1.0_01.count_table',
        dist=rules.prep_dist.output.dist
    output:
        tsv='results/miseq_1.0_01/{method}/miseq_1.0_01.sensspec'
    params:
        outdir='results/miseq_1.0_01/{method}/'
    log:
        'log/miseq_1.0_01/sensspec.method_{method}.miseq_1.0_01.txt'
    shell:
        """
        mothur '#set.logfile(name={log}); set.dir(output={params.outdir});
            sens.spec(list={input.list}, count={input.count_table}, column={input.dist}) '
        """
