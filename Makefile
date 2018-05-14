dag: dag/dag.svg

dag/dag.svg: Snakefile
	snakemake --forceall --rulegraph \
	output/050_deseq/filtered_dds.Rds \
	output/070_clustering/tfs/annotated_clusters_scaled_l2fc.csv \
	output/070_clustering/tfs/annotated_clusters_scaled_l2fc.csv \
	| dot -Tsvg \
	> dag/dag.svg
