set -euo pipefail


## go to the analysis -> kraken2 folder.
module load python/3.9.5

kraken-biom *report.txt --fmt tsv --gzip -o biom_table.tsv

module unload python/3.9.5