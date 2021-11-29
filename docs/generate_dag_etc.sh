# run these commands from repo root to create the DAG, rulegraph and filegraph saved here

# the default rule 'all', plus the three rules 'descriptive', 'missingdata' and 'coinertia_etc'
# (which do not have output files), should be sufficient to re-run the code in its entirety:
snakemake -j 1 -n -F all descriptive missingdata coinertia_etc

# generate DAG etc using:
snakemake -j 1 -n -F all descriptive missingdata coinertia_etc --dag | dot -Tpdf >docs/snakemake_dag.pdf
snakemake -j 1 -n -F all descriptive missingdata coinertia_etc --rulegraph | dot -Tpdf >docs/snakemake_rulegraph.pdf
snakemake -j 1 -n -F all descriptive missingdata coinertia_etc --filegraph | dot -Tpdf >docs/snakemake_filegraph.pdf
