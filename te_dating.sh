#!/bin/bash

DIR=$(cat <<'BABEL_TABLE'
/work/maa146/shane.gossypiumD/
BABEL_TABLE
)
PATH=$PATH:$DIR/bin:$DIR/bin/repeatexplorer/tgicl_linux/bin:/usr/local/R-3.2.3/bin/

R -e 'sessionInfo()'
#R#  version 3.2.3 (2015-12-10)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: CentOS release 6.6 (Final)

# locale:
#  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C
#  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8
#  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8
#  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C
#  [9] LC_ADDRESS=C               LC_TELEPHONE=C
# [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C

# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base


RCMD=$(cat <<'EOF'
args <- commandArgs(trailingOnly = TRUE)

# READ HISTOGRAM
d<-read.table(args[1])

# Fit linear and quadratic model
line<-lm(V2~V1, data=d)
quad<-lm(V2~poly(V1,2), data=d)

l <- quad
category <- 0

# Find best fit using Bayesian Information Criterion and categorize repeat
# cluster based on doi:10.1186/s12864-016-3234-9

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

# Graph histogram and best fit model
png(file.path(dirname(args[1]), "TE_dating.png"))
plot(d)
lines(d$V1, fitted(l))
invisible(dev.off())

# Output data for table
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


# Run mgblast exactly like repeatexplorer on the reads for each cluster, create
# a histogram of % identity, and categorize them with R script above
find $DIR/analysis/clustering/seqClust/clustering/clusters/ -name reads.fas |
    sort |
    while read file; do
        formatdb -p F -i $file
        cat $file |
            parallel --block 100k --recstart '>' --pipe \
                     mgblast -d $file -p85 -W18 -UT -X40 -KT -JF -F '"m D"' -v100000000 -b100000000 -D4 -C 55 -H30 |
            TE_dating_histogram.pl > $(dirname $file)/TE_dating.hist
        Rscript  - <<<"$RCMD"  $(dirname $file)/TE_dating.hist
    done |
    awk 'BEGIN {print "Cluster\ta\tb\tc\tOptimum\tCategory"} 
        {print}' OFS="\t" \
        > $DIR/analysis/clustering/seqClust/clustering/TE_dating.txt
