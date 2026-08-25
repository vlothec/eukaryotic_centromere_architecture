#!/bin/bash

set -e
set -u
set -o pipefail

bin_path_=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
script_dir="$(cd "$(dirname "$0")" && pwd)"

## load lib
source $bin_path_/wfun_lib.sh

filter_comb_="TRUE_FALSE"
min_gap_length_=250
max_gap_length_=100000

while read -r root_name_; do
  printf "\n### $root_name_ ###\n"
  # declare var
  arrays_file_=$(ls data_/${root_name_}*_arrays_expanded.csv 2>/dev/null || true)
  repeats_file_=$(ls data_/${root_name_}*_repeats_filtered.csv 2>/dev/null || true)
  TEanno_file_=$(ls data_/${root_name_}*_edta_filtered.csv.reassigned 2>/dev/null || true)
  intact_file_=$(ls data_/${root_name_}*intact.gff3 2>/dev/null || true)
  karyo_file_=$(ls data_/${root_name_/%.*/}*_karyotype_ext_ 2>/dev/null || true)
  # check they are in place
  [[ -e "$arrays_file_" && -e "$repeats_file_" && -e "$TEanno_file_" && -e "$intact_file_" && -e "$karyo_file_" ]] || { printf "$root_name_: missing files!\n"; continue; }

  filter_comb_="TRUE_FALSE"
  min_gap_length_=250
  max_gap_length_=100000

  # call up the wrapper
  $bin_path_/centromere_dissection_main_ -a $arrays_file_ -r $repeats_file_ -t $TEanno_file_ -i $intact_file_ -f $filter_comb_ -k $karyo_file_ -m $min_gap_length_ -x $max_gap_length_ -n $root_name_ > ${root_name_}_centromere_dissection_main_.log 2>&1 || { printf "centromere_dissection_main_ failed!\n"; continue; }
done < <(ls data_/*_arrays_expanded.csv | awk -F/ '{gsub("_arrays_expanded.csv","",$NF); gsub(".fasta","",$NF); gsub(".fa","",$NF); print $NF}')

