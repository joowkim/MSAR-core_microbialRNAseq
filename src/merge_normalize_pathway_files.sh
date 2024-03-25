set -eou pipefail

module load HUMAnN/3.0

humann_array=($(ls -d1 *out))

mkdir -p humann_merged_results
for i in "${humann_array[@]}" ; do

    sample_name=($(echo $i | sed 's/_humann_out//'))
    humann_join_tables --input $i --output humann_merged_results/$sample_name\_merged_pathway.tsv
    humann_renorm_table --input humann_merged_results/$sample_name\_merged_pathway.tsv --output humann_merged_results/$sample_name\_merged_normalized_pathway.tsv
    humann_split_stratified_table -i humann_merged_results/$sample_name\_merged_normalized_pathway.tsv -o humann_merged_results

done
echo "done!"

module unload HUMAnN/3.0