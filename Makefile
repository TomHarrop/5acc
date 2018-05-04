dag: dag/dag.svg

dag/dag.svg: Snakefile
	snakemake --forceall --rulegraph \
	output/030_mapping/star-pass2/osj/SM_1.ReadsPerGene.out.tab \
	output/010_data/shuffle/shuffed.gtf \
	| dot -Tsvg \
	> dag/dag.svg
