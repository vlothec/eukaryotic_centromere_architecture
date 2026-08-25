#!/bin/bash

set -e
set -u
set -o pipefail

## usage
usage() {
    cat <<EOF

    USAGE
      /path/to/$(basename $0) /path/to/genome.fa /path/to/ann.file.gff3 db app

    DESCRIPTION
      script for running TEsorter after extracting nt sequences of full-llength intact elements

        -db used by TEsorter tool
        -app accepts either all or strict as valid arguments. The first option will consider all intact elements.
         The second, only golden LTRRT (i.e., full length and no nesting structures)

    DEPENDENCIES
      It requires bedtools and TEsorter

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

  ## read ARG
  genome_file=${@:$OPTIND:1}
  if [ ! -f ${genome_file##*\/} ]; then
    # ln -s ${genome_file} .
    cp -p ${genome_file} genome_file.fa
  fi
  genome_file=$(echo ${genome_file##*\/})

  ann_file=${@:$OPTIND+1:1}
  if [ ! -f ${ann_file##*\/} ]; then
    # ln -s ${ann_file} .
    cp -p ${ann_file} ann_file.gff3
  fi
  ann_file=$(echo ${ann_file##*\/})

  db=${@:$OPTIND+2:1}
  if [ -z "${db}" ]; then
    db="plant"
  fi

  app=${@:$OPTIND+3:1}
  if [ -z "${app}" ]; then
    app="all"
  fi

  if [ "${app}" == "all" ]; then
    ## if all annotated TEs
    bedtools_getfasta genome_file.fa <(awk 'BEGIN{OFS="\t"}$3 != "repeat_region" && $3 != "long_terminal_repeat" && $3 != "target_site_duplication"{print $1, $4-1, $5, $1":"$4".."$5"-intact#"$3, $6, $7}' ann_file.gff3) ${ann_file%.*}.fa
    TEsorter ${ann_file%.*}.fa -db rexdb-${db} -p 10 -nolib > ${ann_file%.*}.fa--TEsorter.log 2>&1
  else
    ## If only LTRRT
    ${bin_path}/ann_ltrrt_structure.sh ann_file.gff3 "F"
    bedtools_getfasta genome_file.fa <(awk 'BEGIN{OFS="\t"}$11 == "fl_LTRRT"{print $1, $4-1, $5, $10, $6, $7}' ann_file.gff3.reanno) ${ann_file}.reanno.fa
    TEsorter ${ann_file}.reanno.fa -db rexdb-${db} -p 10 -nolib > ${ann_file}.reanno.fa--TEsorter.log 2>&1
  fi

  test_rm_ith "genome_file.fa"
  test_rm_ith "ann_file.gff3"
