#!/bin/bash

set -e
set -u
set -o pipefail

  bin_path=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

  ## load lib
  source ${bin_path}/wfun_lib.sh

  ## run TEsosrter with reassigned *_edta_filtered_cls.csv
  ## read arg

  genome_file=${@:$OPTIND:1}
  cp -p ${genome_file} genome_file.fa
  genome_file=$(echo ${genome_file##*\/})

  ann_file=${@:$OPTIND+1:1}
  cp -p ${ann_file} ann_file.gff3
  ann_file=$(echo ${ann_file##*\/})

  sample=$(basename ${genome_file})
  echo "## sample: ${sample%%.*}"

  # awk -F, 'BEGIN{OFS = "\t"}{if($4 != "long_terminal_repeat" && $4 ~ /LTR/ && $4 !~ /non_LTR/ && $5 < $6){print $2, $5-1, $6, $2":"$5".."$6"-"$15"#"$4, $7, $8}else if($4 != "long_terminal_repeat" && $4 ~ /LTR/ && $4 !~ /non_LTR/){print $2, $5-1, $6, $2":"$5".."$6"-"$15"#"$4, $7, $8}}' ann_file.gff3 \
    # | awk 'BEGIN{OFS = "\t"}$2 > 0{sub("homology","fragment",$4); sub("structural","intact",$4); print}' \
    # | sort -k1,1 -k2,2n \
    # | awk '$4 ~ /_LTR_|#LTR_/' > ann_file.bed

  awk 'BEGIN{OFS="\t"}{if($0 !~ /^#/ && $3 != "long_terminal_repeat" && $3 != "target_site_duplication" && $5 >  $4){print $0}else if($0 !~ /^#/ && $3 != "long_terminal_repeat" && $3 != "target_site_duplication" && $5 < $4){print $1, $2, $3, $5, $4, $6, $7, $8, $9}}' ann_file.gff3 \
    | awk 'BEGIN{OFS="\t"}{if($0 ~ /homology/){print $1, $4-1, $5, $1":"$4".."$5"-fragment#"$3, $6, $7}else{print $1, $4-1, $5, $1":"$4".."$5"-intact#"$3, $6, $7}}' \
    | awk '$4 ~ /_LTR_|#LTR_/' > ann_file.bed

  printf "\n## TE counts: \n"
  awk '{split($4,arr,"#"); print arr[2]}' ann_file.bed | var_count -

  ## set a min length
  min_length=${@:$OPTIND+3:1}
  if [ -z "${min_length}" ]; then
    min_length=1
  fi

  printf "\n## TEs filtered out after applying min length cutoff ($min_length bp): \n"
  awk -v min_length=${min_length} '$3-$2 < min_length{split($4,arr,"#"); print arr[2]}' ann_file.bed | var_count -
  printf "\n## TEs kept: \n"
  awk -v min_length=${min_length} '$3-$2 >= min_length{split($4,arr,"#"); print arr[2]}' ann_file.bed | var_count -
  awk -v min_length=${min_length} '$3-$2 >= min_length' ann_file.bed > tmp_01
  mv tmp_01 ann_file.bed

  printf "\n\n"
  bedtools_getfasta genome_file.fa ann_file.bed ${ann_file}.fa

  db=${@:$OPTIND+2:1}
  if [ -z "${db}" ]; then
    db="plant"
  fi

  TEsorter ${ann_file}.fa -db rexdb-${db} -p 10 -nolib > ${ann_file}.fa--TEsorter.log 2>&1

