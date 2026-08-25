#!/usr/bin/env bash

set -e
set -u
set -o pipefail

bin_path=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
script_dir="$(cd "$(dirname "$0")" && pwd)"

## load lib
source ${bin_path}/wfun_lib.sh

rescued_file_=${@:$OPTIND:1}
anno_file_=${@:$OPTIND+1:1}

echo "${rescued_file_##*/}"

awk -F "\t" 'BEGIN{OFS="\t"}$0 !~ /^#/{print $3}' ${anno_file_} | var_count - > original_tally_
awk -F "\t" 'BEGIN{OFS="\t"}$0 !~ /^#/{print $3}' ${rescued_file_} | var_count - > reassigned_tally_

awk -F "\t" 'NR==FNR { lookup[$2]=$1; next }{ print $0, ( $2 in lookup ? lookup[$2] : "NA" ) }' original_tally_ reassigned_tally_ \
  | awk 'BEGIN{OFS="\t"}{print $2, $NF, $1, ($1-$NF)}' \
  | sort -k1,1 > merged_tally_

awk 'BEGIN{OFS="\t"; OFMT="%.2f"} NR==FNR { total_original_ += $2; total_rescued_ += $3; total_diff_ += $4; next }{print $0, ($2/total_original_)*100, ($3/total_rescued_)*100, ($4/total_rescued_)*100}' merged_tally_ merged_tally_ > tmp000_;
mv tmp000_ merged_tally_

join -1 1 -2 2 -a 1 -a 2 -e NA -o auto <(sort -k1,1 merged_tally_) <(sort -k2,2 original_tally_) | awk '$2 = $NF' | awk 'NF--' | tr ' ' '\t' > tmp000_;
mv tmp000_ merged_tally_

awk -F "\t" 'NR==FNR { lookup[$3]=$1; next }{ print $0, ( $1 in lookup ? lookup[$1] : "NA" ) }' ${bin_path}/TE_class_list_updated_ merged_tally_ \
  | awk 'BEGIN{OFS="\t"; OFMT="%.2f"}{print $NF, $0}' | awk 'BEGIN{OFS="\t"; OFMT="%.2f"}NF--' | awk 'BEGIN{OFS="\t"; OFMT="%.2f"}{if($1 == "NA"){$1 = $2; print}else{print}}' | sort -k1,1 > tmp000_;
mv tmp000_ merged_tally_

printf "group\tsuperfamily\toriginal\treassigned\tdifference\toriginal_percentage\treassigned_percentage\tdifference_percentage\n" \
  | cat - merged_tally_ > ${rescued_file_##*/}.sum

test_rm_ith original_tally_
test_rm_ith reassigned_tally_
test_rm_ith merged_tally_
