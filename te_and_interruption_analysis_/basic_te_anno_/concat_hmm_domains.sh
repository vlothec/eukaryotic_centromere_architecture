#!/bin/bash

set -e
set -u
set -o pipefail

## usage
usage() {
    cat <<EOF

    USAGE
      /path/to/$(basename $0) /path/to/*.rexdb-plant.cls.tsv /path/to/*.rexdb-plant.dom.faa min_hmm dom /path/to/*.cent.coord

    DESCRIPTION
      script for concatenating hmm functional domains annotated by TEsorter
        -min_hmm: minimum number of hmm-domains for a TE to be considered
        -dom: it takes either all or the specific domain combination you want to concatenate (e.g., GAG|PROT; you must list them
         using a pipe!)
        -To tag the TEs with centromere occurence information, a tab delimited file may be provided, including only chr start end
         and no header

EOF
      exit 0
}

  if [ $# -eq 0 ]; then
    usage
    exit 1
  fi

  bin_path=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

  ## load lib
  source ${bin_path}/wfun_lib.sh

  ## read arg
  export TEsorter_cls_file=${@:$OPTIND:1}
  export TEsorter_aa_file=${@:$OPTIND+1:1}
  export min_hmm=${@:$OPTIND+2:1}
  export dom=${@:$OPTIND+3:1}
  export centromere_ann=${@:$OPTIND+4:1}

  name=$(echo ${TEsorter_cls_file%%.*} | awk -F/ '{print $NF}')

  cp -p ${TEsorter_cls_file} TEsorter_cls_file.tsv
  TEsorter_cls_file="${TEsorter_cls_file##*/}"

  cp -p ${TEsorter_aa_file} TEsorter_aa_file.faa
  TEsorter_aa_file="${TEsorter_aa_file##*/}"

  printf "\n## Sample: ${name}\n"
  printf "\n## min number of HMM: ${min_hmm}\n"

  sed 's/ /:/g' TEsorter_cls_file.tsv > tmp_01
  grep -v 'none$' tmp_01 \
    | paste -d '\t' - <(grep -v 'none$' tmp_01 \
                         | cut -f7 \
                         | awk -F '[:|]' '{for(i=2;i<=NF;i+=2){printf $i" "}; printf "\n" }' \
                         | awk '{ for (i=2; i<=NF; ++i) if ($i != $1) { print "F"; next } print "T" }') \
    | awk '$NF == "T"' | awk 'NF--' \
    | awk 'BEGIN{OFS="\t"}{hmm = split($7,arr,":");{hmm_comb=sprintf("%*s",hmm,"");gsub(/ /,"1",hmm_comb); void_hmm=sprintf("%*s",6-hmm,"");gsub(/ /,"0",void_hmm); $1=$1"/"$3"_"$4"/"hmm_comb""void_hmm; print}}' \
    | awk -v min_hmm=${min_hmm} 'BEGIN{OFS="\t"}{hmm = split($7,arr,":"); if(hmm >= min_hmm){print}}' \
    | awk '$1 ~ $2' > ${TEsorter_cls_file}.5hmm.consistent

  printf "\n## cls summary:\n"
  awk '{print $3"_"$4}' ${TEsorter_cls_file}.5hmm.consistent | var_count -

  printf "\n## Top 5 HMM-dom combinations:\n"
  awk '{print $NF}' ${TEsorter_cls_file}.5hmm.consistent | var_count - | head -n5

  ## if cent.occ is needed
  ## snipet for annotating centromere occurence
  if [ -n "${centromere_ann}" ]; then
    awk 'BEGIN{OFS="\t"}{sub(/\-fra/,"#fra",$1); sub(/\-int/,"#int",$1); print}' ${TEsorter_cls_file}.5hmm.consistent \
      | awk 'BEGIN{OFS="\t"}{split($1,loc,"#|:|\\."); print loc[1], loc[2], loc[4], $0}' \
      | sort -k1,1 \
      | join -1 1 -2 1 -a 1 -a 2 - <(sort -k1,1 ${centromere_ann} | tr ' ' '\t') \
      | tr ' ' '\t' \
      | awk '$NF ~ /[0-9]+/ && NF > 3' \
      | awk 'BEGIN{OFS="\t"}{if($2 >= $(NF-1) && $3 <= $NF){print $0, "IN"}else{print $0, "OUT"}}' \
      | sort -k 13,13 \
      | awk '!a[$4]++' \
      | sort -k 1,1 -k 2,2n \
      | awk 'BEGIN{OFS="\t"}{$4 = $4"/"$NF; print  $4, $5, $6, $7, $8, $9, $10}' \
      | awk 'BEGIN{OFS="\t"}{sub(/#fra/,"-fra",$1); sub(/#int/,"-int",$1); print}' > tmp_01
    mv tmp_01 ${TEsorter_cls_file}.5hmm.consistent

    ## IN/OUT counts
    printf "\n## IN/OUT counts:\n"
    cut -f1 ${TEsorter_cls_file}.5hmm.consistent \
      | awk -F/ '{print $NF}' | var_count -

    printf "\n## IN/OUT counts per cls:\n"
    awk '{sub(/.*\//,"", $1); print $3"_"$4"/"$1}' ${TEsorter_cls_file}.5hmm.consistent | var_count -
  fi

  test_rm_ith "${TEsorter_aa_file}.reanno"
  test_rm_ith "${TEsorter_aa_file}.reanno.concat"
  test_rm_ith "${TEsorter_aa_file}.reanno.concat.${dom/|/_}"
  if [ "$dom" == "all" ]; then
    while read id; do
      grep -A1 ${id%%\#*} TEsorter_aa_file.faa | fasta_to_one - | awk -v id=${id} '{split($1,data,"|"); split(data[2],dom,":"); sub(/.*-/,"",dom[2]); print ">"id"/"dom[2]"\n"$NF}' >> ${TEsorter_aa_file}.reanno
      grep -A1 ${id%%\#*} TEsorter_aa_file.faa | fasta_to_one - | awk 'BEGIN{ORS = ""}{print $NF}' | awk -v id=${id} '{print ">"id"\n"$0}' >> ${TEsorter_aa_file}.reanno.concat.all
    done < <(cut -f1 ${TEsorter_cls_file}.5hmm.consistent)
  else
    while read id; do
      grep -A1 ${id%%\#*} TEsorter_aa_file.faa | fasta_to_one - | awk -v id=${id} '{split($1,data,"|"); split(data[2],dom,":"); sub(/.*-/,"",dom[2]); print ">"id"/"dom[2]"\n"$NF}' >> ${TEsorter_aa_file}.reanno
      grep -A1 ${id%%\#*} TEsorter_aa_file.faa | fasta_to_one - | awk -v dom="${dom}" 'BEGIN{ORS = ""}$1 ~ dom{print $NF}' | awk -v id=${id} '{print ">"id"\n"$0}' >> ${TEsorter_aa_file}.reanno.concat.${dom//|/_}
    done < <(cut -f1 ${TEsorter_cls_file}.5hmm.consistent)
  fi
  test_rm_ith "tmp_01"
  test_rm_ith "TEsorter_aa_file.faa"
  test_rm_ith "TEsorter_cls_file.tsv"
