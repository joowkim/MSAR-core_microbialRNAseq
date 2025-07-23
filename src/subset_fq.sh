set -euo pipefail

output_dir=$2

sample_id=$(echo "$1" | sed 's/[_-]R1\.fastq\.gz//')
new_sample_id=$sample_id\_subset

zcat $sample_id\-R1.fastq.gz | head -n $((100000000 * 4)) | gzip > $output_dir\/$new_sample_id\-R1.fastq.gz &
zcat $sample_id\-R2.fastq.gz | head -n $((100000000 * 4)) | gzip > $output_dir\/$new_sample_id\-R2.fastq.gz &

wait