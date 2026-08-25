## TEsorter_with_EDTA_intact.sh genome.fa ann.file.gff3 db app 
## db should be either plant or metazoa
## app either all or strict (the former considering all intact elements the latter golden LTRRT)

## TEsorter_diagnostic.sh *.rexdb-plant.cls.tsv all
## either all or the lineages you want to look at

## concat_hmm_domains.sh *.rexdb-plant.cls.tsv *.rexdb-plant.dom.faa min_hmm dom *.cent.coord
## min_hmm, minimum number of hmm domains for a TE be considered
## dom might take either all or the specific domain combination you want to concatenate e.g., GAG|PROT

## subset_concat_hmm_domains.sh *.rexdb-plant.cls.tsv.5hmm.consistent "GAG|CRM:PROT|CRM:RT|CRM:RH|CRM:INT|CRM:CHDCR|CRM,GAG|CRM:PROT|CRM:RT|CRM:RH|CRM:INT|CRM” 
## second argument must include a comma separated value with target domain combinations

## mafft_plus_fasttree_alone.sh *.rexdb-plant.dom.faa.reanno.concat.all cls alg maxiterate
## specifying cls is necessary to pull out an appropriate root sequence (e.g.,Gypsy_Athila)
## mafft alignment: retree or globalpair
## maxiterate is an argument for mafft run (e.g., 100)


## TEsorter_with_filtered_rescued_files.sh genome.fa filtered.rescue.ann.file.csv db min_length 
