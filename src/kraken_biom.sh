set -euo pipefail


## go to the analysis -> kraken2 folder.
source ~/beegfs/python_env/kraken2/bin/activate

kraken-biom *report.txt --fmt tsv --gzip -o biom_table.tsv

deactivate