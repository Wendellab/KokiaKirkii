#!/bin/bash

DIR=$(cat <<'BABEL_TABLE'
/scratch/maa146/shane.gossypiumD/
BABEL_TABLE
)
PATH=$PATH:$DIR/bin:$DIR/bin/repeatexplorer/tgicl_linux/bin:/usr/local/R-3.2.3/bin/

RCMD=$(cat <<'EOF'
args <- commandArgs(trailingOnly = TRUE)

d<-read.table(args[1])
line<-lm(V2~V1, data=d)
quad<-lm(V2~poly(V1,2), data=d)

l <- quad
category <- 0

if(BIC(line) < BIC(quad)){
    l <- line
    if(l$coefficients[2] > 0.001){
        category <- 1
    } else if (l$coefficients[2] < 0.001 && l$coefficients[2] > -0.001)  {
        category <- 2
    } else if (l$coefficients[2] < -0.001) {
        category <- 3
    }
} else {
    opti <- d[which.max(fitted(l)),1]
    
    if(l$coefficients[3] > 0){
        category <- 4
        if(l$coefficients[2]<0){
            category <- "4*"
        }
    } else  {
        if(opti < 99) {
            category <- 5
        } else {
            category <- 6
        }

    }

}

png(file.path(dirname(args[1]), "TE_dating.png"))
plot(d)
lines(d$V1, fitted(l))
invisible(dev.off())

cluster <- sub(".*dir_", "", dirname(args[1]))
cat(sprintf("%s\t%f\t%f\t%f\t%d\t%s\n",
            cluster,
            l$coefficients[1],
            l$coefficients[2],
            l$coefficients[3],
            d[which.max(fitted(l)),1],
            category))
EOF
)


find $DIR/analysis/clustering/seqClust/clustering/clusters/ -name reads.fas |
    sort |
    while read file; do
        # formatdb -p F -i $file
        # cat $file |
        #     parallel --block 100k --recstart '>' --pipe \
        #              mgblast -d $file -p85 -W18 -UT -X40 -KT -JF -F '"m D"' -v100000000 -b100000000 -D4 -C 55 -H30 |
        #     TE_dating_histogram.pl > $(dirname $file)/TE_dating.hist
        Rscript  - <<<"$RCMD"  $(dirname $file)/TE_dating.hist
    done |
    awk 'BEGIN {print "Cluster\ta\tb\tc\tOptimum\tCategory"} 
        {print}' OFS="\t" \
        > $DIR/analysis/clustering/seqClust/clustering/TE_dating.txt
