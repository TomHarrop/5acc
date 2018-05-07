dag: dag/dag.svg

dag/dag.svg: Snakefile
	snakemake --forceall --rulegraph \
	output/040_background-counts/all_counts.Rds \
	output/060_tpm/tpm.Rds \
	| dot -Tsvg \
	> dag/dag.svg
