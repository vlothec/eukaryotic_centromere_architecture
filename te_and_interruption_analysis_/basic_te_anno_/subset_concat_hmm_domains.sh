#!/bin/bash

set -e
set -u
set -o pipefail

## usage
usage() {
    cat <<EOF

    USAGE
      /path/to/$(basename $0) /path/to/*.rexdb-plant.cls.tsv.5hmm.consistent dom_comb

    DESCRIPTION
      It generates a list including only the TEs meeting any of the domain combinations provided in the second ARG
        -dom_comb: a comma-separated string quoted list
         (e.g., "GAG|CRM:PROT|CRM:RT|CRM:RH|CRM:INT|CRM:CHDCR|CRM,GAG|CRM:PROT|CRM:RT|CRM:RH|CRM:INT|CRM”)
         Exactly as they appear in column 7 of *.rexdb-plant.cls.tsv file.

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

  export TEsorter_cls_file=${@:$OPTIND:1}
  export dom_list=${@:$OPTIND+1:1}
  dom_list=( $(echo $dom_list | tr ',' ' ') )

  if [ ! -f ${TEsorter_cls_file}.concat.subset.list.1 ]; then
    id=1
  else
    id=$(ls ${TEsorter_cls_file}.concat.subset.list.* | awk -F "." '{print $NF}' | sort -k1,1nr | awk 'NR==1{print $0+1}')
  fi

  test_rm_ith concat.subset.list
  for dom in "${dom_list[@]}"; do
    awk -v dom="${dom}" '$7 == dom{print $1}' ${TEsorter_cls_file} >> ${TEsorter_cls_file}.concat.subset.list.${id}
  done


