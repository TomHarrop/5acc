dag: dag/dag.svg

dag/dag.svg: Snakefile
	snakemake --forceall --rulegraph \
	output/100_figures/Figure_1.pdf \
	output/100_figures/Figure_2.pdf \
	output/100_figures/Figure_3.pdf \
	output/100_figures/Figure_4.pdf \
	output/100_figures/Figure_5.pdf \
	output/100_figures/Figure_6.pdf \
	output/100_figures/Figure_S2.pdf \
	output/100_figures/Figure_S5.pdf \
	output/100_figures/Figure_S6.pdf \
	output/100_figures/Figure_S8.pdf \
	output/100_figures/Figure_S9.pdf \
	output/110_tables/Table_S4.csv \
	output/110_tables/Table_S6.csv \
	output/110_tables/Table_S8.csv \
	output/110_tables/Table_S9.csv \
	| dot -Tsvg \
	> dag/dag.svg
