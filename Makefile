dag: dag/dag.svg

dag/dag.svg: Snakefile
	snakemake --forceall --rulegraph \
	output/050_deseq/filtered_dds.Rds \
	| dot -Tsvg \
	> dag/dag.svg
