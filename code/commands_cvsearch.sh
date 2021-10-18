bash code/process_mice.sh data/miseq/miseq.files
mothur "#deunique.seqs(fasta=data/miseq/miseq.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta,  count=data/miseq/miseq.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.pick.pick.count_table)"
rm data/miseq/miseq.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.pick.pick.redundant.groups
sed "s/_/-/g" < data/miseq/miseq.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.redundant.fasta > data/miseq/miseq.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.redundant.fix.fasta
R -e "source('code/generate_samples.R'); generate_indiv_samples('data/miseq/miseq.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.redundant.fix.fasta', 'data/miseq/miseq', 1.0, '01')"
mothur "#unique.seqs(fasta=data/miseq/miseq_1.0_01.fasta)"
mothur "#dist.seqs(fasta=data/miseq/miseq_1.0_01.unique.fasta, processors=8, cutoff=0.20)"
mothur "#degap.seqs(fasta=data/miseq/miseq_1.0_01.fasta)"
bash code/run_cvsearch.sh data/miseq/miseq_1.0_01.ng.fasta
mothur "#sens.spec(column=data/miseq/miseq_1.0_01.unique.dist, list=data/miseq/miseq_1.0_01.cvsearch.list, name=data/miseq/miseq_1.0_01.names, label=userLabel, cutoff=0.03, outputdir=data/miseq)"
