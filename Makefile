dag: dag/dag.svg

dag/dag.svg: Snakefile
	snakemake --forceall --rulegraph \
	output/060_tpm/tpm_with_calls.Rds \
	| dot -Tsvg \
	> dag/dag.svg
