dag: dag/dag.svg

dag/dag.svg: Snakefile
	snakemake --forceall --rulegraph \
	output/050_deseq/wald_tests/sig/domestication.csv \
	output/050_deseq/tfs/sig/domestication.csv \
	output/070_clustering/tfs/annotated_clusters_scaled_l2fc.csv \
	output/070_clustering/all/annotated_clusters_scaled_l2fc.csv \
	| dot -Tsvg \
	> dag/dag.svg
