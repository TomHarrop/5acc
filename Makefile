dag: text/dag.svg

text/dag.svg: Snakefile
	snakemake --forceall --dag \
	output/030_mapping/star-pass2/osj/SM_1.ReadsPerGene.out.tab \
	| dot -Tsvg \
	> text/dag.svg