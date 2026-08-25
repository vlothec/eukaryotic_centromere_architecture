#!/bin/bash

set -e
set -u
set -o pipefail

  bin_path=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

  ## load lib
  source ${bin_path}/wfun_lib.sh

  # set var

  gap_coordinates_file_=${@:$OPTIND:1}
  TEanno_file_=${@:$OPTIND+1:1}
  root_name_=${@:$OPTIND+2:1}
  filter_comb_=${@:$OPTIND+3+1}

  # sort _gap_coordinates_file_.tsv input
  # -1 to make it 0 based
  printf "$gap_coordinates_file_\n"
  awk 'NR>1' ${gap_coordinates_file_} \
    | awk 'BEGIN{OFS="\t"}{if($NF !~ /_/){gap_id_ = NR"_"$NF; print $1, $4, $5, gap_id_}else{print}}' \
    | awk 'BEGIN{OFS="\t"}{$2 = $2-1; print}' \
    | awk 'BEGIN{OFS="\t"}{if($2 < 0){$2 = 0; print}else{print}}' \
    | sort -k1,1 -k2,2n -k3,3n > ${root_name_}_sorted_gap_coordinates_file_

  # sort edta_filtered.csv.reassigned input
  printf "$TEanno_file_\n"
  awk -F, '$0 !~ /repeat_region/ && $0 !~ /repeat_region*structural/ && $0 !~ /long_terminal_repeat|target_site_duplication/ && $0 !~ /Sequence_ontology/' ${TEanno_file_} \
    | awk -F, 'BEGIN{OFS="\t"}{print $2, $5, $6, $4, $7, $8}' \
    | awk 'BEGIN{OFS="\t"}{if($2 < 0){$2 = 0; print}else{print}}' \
    | awk 'BEGIN{OFS="\t"}{$2=sprintf("%.0f",$2); $3=sprintf("%.0f",$3); print}' \
    | sort -k1,1 -k2,2n -k3,3n > ${root_name_}_sorted_TEanno_file_
  
  wc -l ${root_name_}_sorted_TEanno_file_
  printf "${root_name_}_sorted_TEanno_file_\n"
  
  # intersect
  # max_coord_=$(awk '{print $3}' ${root_name_}_sorted_TEanno_file_ | sort -k1,1nr - | head -n1 -)
  max_coord_=1
  printf "max coord found: $max_coord_\n"
  
  if (( max_coord_ < 512000000 )); then
    printf "bedtools will be used\n"
    bedtools intersect -a ${root_name_}_sorted_gap_coordinates_file_ -b ${root_name_}_sorted_TEanno_file_ -wao \
      | awk -v filter_comb_=$filter_comb_ 'BEGIN{OFS="\t"}{print $0, filter_comb_}' > ${root_name_}_TEanno_isec_gaps_
  else
    printf "bedmaps will be used\n"
    bedmap --echo --echo-map --bases-uniq --delim '\t' ${root_name_}_sorted_gap_coordinates_file_ ${root_name_}_sorted_TEanno_file_ \
      | awk 'BEGIN{OFS="\t"}$NF == 0{print $1, $2, $3, $4, ".", "-1", "-1", ".", "-1", ".", $NF}' > ${root_name_}_TEanno_isec_gaps_
    bedmap --echo --echo-map --bases --delim '\t' ${root_name_}_sorted_TEanno_file_ ${root_name_}_sorted_gap_coordinates_file_ \
      | awk 'BEGIN{OFS="\t"}{if(NF == 11 && $NF > 0){print $7, $8, $9, $10, $1, $2, $3, $4, $5, $6, $NF}}' >> ${root_name_}_TEanno_isec_gaps_
    sort -k1,1 -k2,2n ${root_name_}_TEanno_isec_gaps_ \
       | awk -v filter_comb_=$filter_comb_ 'BEGIN{OFS="\t"}{print $0, filter_comb_}' > ${root_name_}_TEanno_isec_gaps_tmp000
    mv ${root_name_}_TEanno_isec_gaps_tmp000 ${root_name_}_TEanno_isec_gaps_
  fi

  rm -I "${root_name_}_sorted_TEanno_file_"
  rm -I "${root_name_}_sorted_gap_coordinates_file_"
