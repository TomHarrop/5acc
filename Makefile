dag: dag/dag.svg

dag/dag.svg: Snakefile
	snakemake --forceall --rulegraph \
	output/040_background-counts/osj/SM_1.htseq-count \
	output/030_mapping/stats/star_logs.csv \
	output/050_deseq/dds.Rds \
	| dot -Tsvg \
	> dag/dag.svg
