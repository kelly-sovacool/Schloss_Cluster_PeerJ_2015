perc_identity = 0.97  # to match mothur's 0.3 dissimilarity threshold
min_seq_length = 30 # from Pat's vsearch script
max_accepts = 16
max_rejects = 64
# the default value of wordlength is already 8 but I'm paranoid
word_length = 8

datasets = ["miseq_PDS", "mouse_KLS"]

rule targets:
    input:
        expand("results/{dataset}/de_novo/{dataset}.sensspec",
                dataset = datasets),
        'results/seq_counts.txt'

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

rule copy_pds_files:
    input:
        count_table='data/miseq/miseq.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.pick.pick.count_table',
        fasta='data/miseq/miseq.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta'
    output:
        count_table="data/miseq_PDS/miseq_PDS.count_table",
        fasta="data/miseq_PDS/miseq_PDS.fasta"
    shell:
        """
        cp {input.count_table} {output.count_table}
        cp {input.fasta} {output.fasta}
        """

rule calc_dists_dataset:
    input:
        fasta="data/{dataset}/{dataset}.fasta"
    output:
        column="data/{dataset}/{dataset}.dist"
    params:
        outdir="data/{dataset}/"
    log:
        "log/{dataset}/calc_dists.log"
    resources:
        procs=16
    shell:
        """
        mothur '#set.logfile(name={log}); set.dir(output={params.outdir});
            dist.seqs(fasta={input.fasta}, cutoff=0.03, processors={resources.procs}) '
        """

# replace dots/hyphens in fasta headers (1st & 2nd occurrence),
# but not in the distance value (3rd occurrence).
rule prep_dist:
    input:
        dist="data/{dataset}/{dataset}.dist"
    output:
        dist='data/{dataset}/{dataset}.ng.dist'
    shell:
        """
        cat {input.dist} |  sed 's/[\.-]/_/' | sed 's/[\.-]/_/' > {output.dist}
        """

# vsearch doesn't support dots or hyphens in sequence names.
rule prep_count_table:
    input:
        count_table="data/{dataset}/{dataset}.count_table"
    output:
        count_table='data/{dataset}/{dataset}.ng.count_table'
    shell:
        """
        cat {input.count_table} |  sed 's/[\.-]/_/g' > {output.count_table}
        """

rule remove_gaps_query:
    input:
        fna="data/{dataset}/{dataset}.fasta"
    output:
        degap=temp("data/{dataset}/{dataset}.tmp.fasta"),
        fna="data/{dataset}/{dataset}.ng.fasta"
    log:
        'log/{dataset}/remove_gaps_query.log'
    params:
        outdir="data/{dataset}/",
        prefix='{dataset}.tmp'
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
        fna="data/{dataset}/{dataset}.ng.fasta"
    output:
        fna="data/{dataset}/{dataset}.ng.sorted.fasta",
        uc="data/{dataset}/{dataset}.ng.sorted.uc"
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
        uc='results/{dataset}/de_novo/{dataset}.uc'
    benchmark:
        'benchmarks/{dataset}/vsearch.method_de_novo.{dataset}.txt'
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
        code='code/uc_to_list_KLS.R',
        sorted=rules.vsearch_sort.output.uc,
        clustered='results/{dataset}/{method}/{dataset}.uc'
    output:
        list='results/{dataset}/{method}/{dataset}.list'
    script:
        'code/uc_to_list_KLS.R'


rule sensspec_vsearch:
    input:
        list=rules.uc_to_list.output.list,
        count_table=rules.prep_count_table.output.count_table,
        dist=rules.prep_dist.output.dist
    output:
        tsv='results/{dataset}/{method}/{dataset}.sensspec'
    params:
        outdir='results/{dataset}/{method}/'
    log:
        'log/{dataset}/sensspec.method_{method}.{dataset}.txt'
    shell:
        """
        mothur '#set.logfile(name={log}); set.dir(output={params.outdir});
            sens.spec(list={input.list}, count={input.count_table}, column={input.dist}) '
        """

rule count_seqs:
    input:
        expand('data/{dataset}/{dataset}.fasta', dataset = datasets)
    output:
        txt='results/seq_counts.txt'
    shell:
        """
        for f in {input}; do
            nseqs=$(grep '>' $f | wc -l)
            echo "${{f}}\t${{nseqs}}" >> {output.txt}
        done
        """

rule uc_to_list_MISEQ1:
    input:
        code='code/uc_to_list_KLS.R',
        sorted='data/miseq/miseq_1.0_01.vdgc.sorted.uc',
        clustered='data/miseq/miseq_1.0_01.vdgc.clustered.uc'
    output:
        list='results/miseq_1.0_01/de_novo/miseq_1.0_01.list'
    script:
        'code/uc_to_list_KLS.R'

rule sensspec_vsearch_MISEQ1:
    input:
        list=rules.uc_to_list_MISEQ1.output.list,
        count_table="data/miseq_PDS/miseq_PDS.count_table",
        dist='data/miseq/miseq_1.0_01.unique.dist'
    output:
        tsv='results/miseq_1.0_01/de_novo/miseq_1.0_01.sensspec'
    params:
        outdir='rresults/miseq_1.0_01/de_novo/'
    log:
        'log/miseq_1.0_01/sensspec.method_de_novo.miseq_1.0_01.txt'
    shell:
        """
        mothur '#set.logfile(name={log}); set.dir(output={params.outdir});
            sens.spec(list={input.list}, count={input.count_table}, column={input.dist}) '
        """