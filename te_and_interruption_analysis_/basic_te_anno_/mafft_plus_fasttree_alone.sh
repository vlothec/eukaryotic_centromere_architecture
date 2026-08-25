#!/bin/bash

set -e
set -u
set -o pipefail

## usage
usage() {
    cat <<EOF

    USAGE
      /path/to/$(basename $0) /path/to/*.rexdb-plant.dom.faa.reanno.concat.all cls alg maxiterate

    DESCRIPTION
      It runs mafft and fasttree
      specifying cls is necessary to pull out an appropriate root sequence (e.g.,Gypsy_Athila)
      mafft alignment: retree or globalpair
      maxiterate is an argument for mafft run (e.g., 100)

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
  export concat_file=${@:$OPTIND:1}
  export cls=${@:$OPTIND+1:1}
  export alg=${@:$OPTIND+2:1}
  export maxiterate=${@:$OPTIND+3:1}

      # append root sequence
      if ! grep -P -q "/Saccharomyces_cerevisiae" ${concat_file}; then
        grep -A1 ${cls%_*} ${bin_path}/sce_outliers.concat_dom.faa | grep -v '\-\-' >> ${concat_file}
      fi

      seq_file="${concat_file}"
      if [ "${alg}" == "retree" ]; then
        mafft --retree 2 --maxiterate ${maxiterate} --thread 10 --memsave ${seq_file} > ${seq_file}.retree2.max${maxiterate}.mafft 2> ${seq_file}.retree2.max${maxiterate}.mafft.log
      fi

      if [ "${alg}" == "globalpair" ]; then
        mafft --globalpair --maxiterate ${maxiterate} --thread 10 ${seq_file} > ${seq_file}.globalpair.max${maxiterate}.mafft 2> ${seq_file}.globalpair.max${maxiterate}.mafft.log
      fi

      sed -i 's/:\|\./\&/g' ${seq_file}.retree2.max${maxiterate}.mafft
      FastTree ${seq_file}.retree2.max${maxiterate}.mafft > ${seq_file}.retree2.max${maxiterate}.mafft.fasttree 2> ${seq_file}.retree2.max${maxiterate}.mafft.fasttree.log

