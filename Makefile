dag: text/dag.pdf text/dag.svg

text/dag.pdf: Snakefile
	snakemake --forceall --dag \
	output/030_mapping/star-pass2/osj/SM_1.ReadsPerGene.out.tab \
	| dot -Tpdf \
	> text/dag.pdf

text/dag.svg: Snakefile
	snakemake --forceall --dag \
	output/030_mapping/star-pass2/osj/SM_1.ReadsPerGene.out.tab \
	| dot -Tsvg \
	> text/dag.svg