dag: dag/dag.svg

dag/dag.svg: Snakefile
	snakemake --forceall --rulegraph \
	output/050_deseq/filtered_dds.Rds \
	output/050_deseq/sig/stage.tab \
	| dot -Tsvg \
	> dag/dag.svg
