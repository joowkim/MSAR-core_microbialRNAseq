#set -euo pipefail


for i in *.profiled_metagenome.txt
do
  sample_name=($(echo $i | sed 's/\.profiled_metagenome\.txt//'))
  cat $i | grep s__ | grep -v t__ | cut -f1,5 | sed "s/.*|//g" > $sample_name.species.txt
  cat $i | grep g__ | grep -v t__ | grep -v s__ | cut -f1,5 | sed "s/.*|//g" > $sample_name.genus.txt

done