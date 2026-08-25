#!/bin/bash

set -e
set -u
set -o pipefail


  ## load functions
  ## get sequences
  bedtools_getfasta () {
    if [ -f ${1}.fai ];then
      rm -I ${1}.fai
    fi
    bedtools getfasta -name -s -fi ${1} -bed ${2} -fo ${3}
    sed -i 's/(.)//g' ${3}
  }

  ## parse fasta files one row/chr
  fasta_to_one () {
    awk 'BEGIN{OFS="\t"}{if($0 ~ ">"){ORS = FS; print "\n"$0}else if($0 !~ ">"){ORS = "";print $0}else{print "\n"}}' ${1} \
      | awk 'NR>1{print $0}'
  }

  ## parse fasta and get seq length
  fasta_seq_length () {
    awk 'BEGIN{OFS="\t"}{if($0 ~ ">"){ORS = FS; print "\n"$0}else if($0 !~ ">"){ORS = "";print $0}else{print "\n"}}' ${1} \
      | awk 'NR>1{print $0}' \
      | awk 'BEGIN{OFS="\t"}{print $1, length($2)}'
  }

  ## remove files in a controlled way
  test_rm_ith () {
    counter=0
    for file in ${1};do
      if [ -f ${file} ];then
        counter=$[$counter +1]
      fi
    done

    ## test how many files will be treated
    if [[ $counter -lt 20 ]];then
      printf "$counter files will be rm:\n" >> test_rm_ith.log
    else
      printf "Too many args. Aborted!" >> test_rm_ith.log
      exit 1
      ## return 1 might be used too as long as you capture it;i.e., (test_rm_ith) || exit $?
    fi

    for file in ${1};do
      if [ -f ${file} ];then
        printf "file $file rm\n" >> test_rm_ith.log
        rm -I ${file}
      fi
    done
    printf " \n" >> test_rm_ith.log
  }

  var_count () {
    sort ${1} | uniq -c | awk '{$1=$1};1' | sort -k 1,1nr | awk 'BEGIN{OFS="\t"}NF == 2{print $1, $2}'
  }

