dag: dag/dag.svg

dag/dag.svg: Snakefile
	snakemake --forceall --rulegraph \
	output/100_figures/Figure_1.pdf \
	output/100_figures/Figure_S2.pdf \
	| dot -Tsvg \
	> dag/dag.svg
