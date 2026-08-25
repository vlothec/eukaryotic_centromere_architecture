#!/bin/bash

set -e
set -u
set -o pipefail


  ## gff3 input file is scanned to identify and flag fl_LTR, lLTR and rLTR

  var_count () {
    sort ${1} | uniq -c | awk '{$1=$1};1' | sort -k 1,1nr | awk 'BEGIN{OFS="\t"}NF == 2{print $1, $2}'
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
    if [[ $counter -lt 10 ]];then
      printf "\nCLEAN-UP: $counter files will be rm\n"
    else
      printf "Too many args. Aborted!"
      exit 1
      ## return 1 might be used too as long as you capture it;i.e., (test_rm_ith) || exit $?
    fi

    for file in ${1};do
      if [ -f ${file} ];then
        rm -I ${file}
      fi
    done
  }

  ann_file=${@:$OPTIND:1}

  cp -p ${ann_file} ann_file.gff3
  ann_file=$(echo ${ann_file##*\/})

  TEsorter=${@:$OPTIND+1:1}

  awk '$3 !~ /repeat_region|target_site_duplication/' ann_file.gff3 \
    | awk 'BEGIN{OFS="\t"}$0 !~ /^ /{print $0, $1"_"$4, $1"_"$5, sqrt(($5-$4)^2)}' \
    | sort -k1,1 -k4,4n -k12,12nr \
    | awk 'BEGIN{OFS="\t"}NF--' > sample.intact.gff3

  ## fl_LTRRT
  awk 'NR==FNR && !/^#/{a[$10]++;b[$11]++;next}a[$10] > 1 && b[$11] > 1 && a[$10] < 3 && b[$11] < 3' sample.intact.gff3 sample.intact.gff3 \
    | awk '{print $0"\t"$1":"$4".."$5"#"$3"\tfl_LTRRT";}' \
    | tr ' ' '\t' \
    | cut -f1-9,12,13 > tmp_01

  ## lLTR
  awk 'NR==FNR && !/^#/{a[$10]++;b[$11]++;next}a[$10] > 1 && a[$10] < 3' sample.intact.gff3 sample.intact.gff3 \
    | awk 'NR==FNR && !/^#/{a[$10]++;b[$11]++;next}b[$11] < 2' - <(awk 'NR==FNR && !/^#/{a[$10]++;b[$11]++;next}a[$10] > 1 && a[$10] < 3' sample.intact.gff3 sample.intact.gff3) \
    | awk 'start == $4 || end == $5{print $0"\t"id}{start = $4; end = $5; id = $1":"$4".."$5"#"$3"\tlLTR";}' \
    | tr ' ' '\t' \
    | cut -f1-9,12,13 >> tmp_01

  ## rLTR
  awk 'NR==FNR && !/^#/{a[$10]++;b[$11]++;next}b[$11] > 1 && b[$11] < 3' sample.intact.gff3 sample.intact.gff3 \
    | awk 'NR==FNR && !/^#/{a[$10]++;b[$11]++;next}a[$10] < 2' - <(awk 'NR==FNR && !/^#/{a[$10]++;b[$11]++;next}b[$11] > 1 && b[$11] < 3' sample.intact.gff3 sample.intact.gff3) \
    | awk 'start == $4 || end == $5{print $0"\t"id}{start = $4; end = $5; id = $1":"$4".."$5"#"$3"\trLTR";}' \
    | tr ' ' '\t' \
    | cut -f1-9,12,13 >> tmp_01
  sort -k 4,4n -k 10,10 -k 11,11 tmp_01 > sample.intact.gff3

  ## remove complex structures; i.e., nested elements
  grep -v -f <(cut -f10 sample.intact.gff3 | var_count - | awk '!/^3/{split($2,arr,"#"); print arr[1]}') sample.intact.gff3 \
    | awk '$10 ~ $4 || $10 ~ $5' > tmp_01
  mv tmp_01 sample.intact.gff3


  ## keep hits with consistent domain annotation
  if [ ${TEsorter} == "T" ]; then
    name="${ann_file##*/}"
    test_rm_ith "sample.intact.fa.rexdb-plant.cls.tsv"
    ln -s ${name%.*}.fa.rexdb-plant.cls.tsv sample.intact.fa.rexdb-plant.cls.tsv
    grep -v 'none$' sample.intact.fa.rexdb-plant.cls.tsv \
      | paste -d '\t' - <(grep -v 'none$' sample.intact.fa.rexdb-plant.cls.tsv \
                           | cut -f7 \
                           | awk -F '[ |]' '{for(i=2;i<=NF;i+=2){printf $i" "}; printf "\n" }' \
                           | awk '{ for (i=2; i<=NF; ++i) if ($i != $1) { print "F"; next } print "T" }') \
      | awk -F '[\t]' 'BEGIN{OFS="\t"}{if($8 == "T"){print $1,$3"_"$4}}' \
      | sed 's/:\|\.\./_/g' \
      | awk 'BEGIN{OFS="\t"}{split($1,coord,"#"); $1 = coord[1]; print}' > sample.cls

    ## add cls data into sample.intact.gff3 file
    ## unclassified LTRRTs will be assigned a generic NA
    sort -k 1,1 sample.cls \
      | join -1 1 -2 12 -a 1 -a 2 - <(awk 'BEGIN{OFS="\t"}{split($10,coord,"#"); print $0, coord[1]}' sample.intact.gff3 | sort -k12,12) \
      | awk 'BEGIN{OFS="\t"}{if(NF == 13){print $3, $4, $5, $6, $7, $8, $9, $10, $11, $12"/"$2, $13}else if(NF == 12){print $2, $3, $4, $5, $6, $7, $8, $9, $10, $11"/NA", $12}}' > tmp_01
    mv tmp_01 sample.intact.gff3
  fi

  mv sample.intact.gff3 ${ann_file##*/}.reanno

  test_rm_ith "sample.intact.fa.rexdb-plant.cls.tsv"
  test_rm_ith "sample.cls"
  test_rm_ith "ann_file.gff3"
