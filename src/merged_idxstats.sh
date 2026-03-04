# 1. Get gene names from first file
first=$(ls DRS_*.idxstat | head -n 1)
cut -f1 "$first" > genes.txt

# 2. Extract mapped reads (3rd column) from all samples
for f in DRS_*.idxstat; do
    cut -f3 "$f" > "$f.count"
done

# 3. Create header
echo -e "Gene\t$(ls DRS_*.idxstat | sed 's/.idxstat//g' | tr '\n' '\t')" > header.txt

# 4. Merge everything
paste genes.txt DRS_*.idxstat.count > body.txt
cat header.txt body.txt > gene_count_matrix.txt