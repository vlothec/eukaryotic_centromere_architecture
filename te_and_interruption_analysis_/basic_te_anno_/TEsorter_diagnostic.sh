#!/bin/bash

set -e
set -u
set -o pipefail

## usage
usage() {
    cat <<EOF

    USAGE
      /path/to/$(basename $0) /path/to/*.rexdb-plant.cls.tsv clade

    DESCRIPTION
      script for summarising TEsorter hmm-domain combinations
        -clade: either all or specific one

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
  export cls=${@:$OPTIND+1:1}

  name=$(echo ${TEsorter_cls_file%%.*} | awk -F/ '{print $NF}')

  cp -p ${TEsorter_cls_file} TEsorter_cls_file.tsv
  TEsorter_cls_file="${TEsorter_cls_file##*/}"

  if [ "$cls" == "all" ]; then

    if [ -f "${TEsorter_cls_file}.all.sum" ];then
      test_rm_ith "${TEsorter_cls_file}.all.sum"
    fi

    while read cls; do
      printf "\n## Clade: ${cls}\n"
      awk -v cls="${cls}" 'tolower($4) == tolower(cls){print $2"_"$3"_"$4}' TEsorter_cls_file.tsv | var_count -
      printf "\n## T0P10 combination domain:\n"
      sed 's/ /:/g' TEsorter_cls_file.tsv | awk -v cls="${cls}" 'tolower($4) == tolower(cls) && $1 ~ $2{print $7}' | var_count - | head -n10
      sed 's/ /:/g' TEsorter_cls_file.tsv \
        | awk -v cls="${cls}" 'tolower($4) == tolower(cls) && $1 ~ $2{print $7}' \
        | var_count - \
        | awk -v file="${TEsorter_cls_file}" -v cls="${cls}" 'BEGIN{OFS="\t"}{n = split($2,dom,":"); print file, cls, $1, n, $2}' >> ${TEsorter_cls_file}.all.sum
      printf "\n\n"
    done < <(awk '{print $4}' TEsorter_cls_file.tsv | sort | uniq)
  else
    printf "\n## Clade: ${cls}\n"
    awk -v cls="${cls}" 'tolower($4) == tolower(cls){print $2"_"$3"_"$4}' TEsorter_cls_file.tsv | var_count -
    printf "\n## T0P10 combination domain:\n"
    sed 's/ /:/g' TEsorter_cls_file.tsv | awk -v cls="${cls}" 'tolower($4) == tolower(cls) && $1 ~ $2{print $7}' | var_count - | head -n10
    sed 's/ /:/g' TEsorter_cls_file.tsv \
      | awk -v cls="${cls}" 'tolower($4) == tolower(cls) && $1 ~ $2{print $7}' \
      | var_count - \
      | awk -v file="${TEsorter_cls_file}" -v cls="${cls}" 'BEGIN{OFS="\t"}{n = split($2,dom,":"); print file, cls, $1, n, $2}' > ${TEsorter_cls_file}.${cls}.sum
    printf "\n\n"
  fi

  test_rm_ith "TEsorter_cls_file.tsv"
