dag: dag/dag.svg

dag/dag.svg: Snakefile
	snakemake --forceall --rulegraph \
	output/050_deseq/wald_tests/expr_genes/sig/domestication.csv \
	output/050_deseq/wald_tests/tfs/sig/domestication.csv \
	output/070_clustering/tfs/annotated_clusters_scaled_l2fc.csv \
	output/080_phenotype/cali_pca.Rds \
    output/080_phenotype/mtp_cluster_correlation.csv \
	| dot -Tsvg \
	> dag/dag.svg
