#!/bin/env Rscript

options(warn=1)

# get script name
all.args = commandArgs(F)
fn.arg = "--file="
script.name = sub(fn.arg, "", all.args[grep(fn.arg, all.args)])

args = commandArgs(T)
if (length(args) == 0) {
  cat(sprintf("usage: %s <contig-bin table> <contigs> <contig field> <ofn>\n", script.name))
  q(status=1)
}

ifn.cb = args[1]
ifn.contigs = args[2]
contig.field = args[3]
ofn = args[4]

cb = read.delim(ifn.cb, header=F)
colnames(cb) = c("contig", "bin")

contigs = read.delim(ifn.contigs)
contigs$contig = contigs[,contig.field]
contigs$length = contigs$end - contigs$start

field.count = function(x, field="gene")
{
    tt = table(x[,field])
    result = data.frame(x=names(tt), count=as.vector(tt))
    names(result)[1] = field
    result[order(result$count, decreasing=T),]
}

bins = sort(unique(cb$bin))
fc = field.count(cb, "contig")
if (!all(fc$count == 1)) {
    cc = fc$contig[fc$count > 1]
    stop(sprintf("contig %s appears more than once in table: %s", cc[1], ifn.cb))
}
cb$length = contigs$length[match(cb$contig,contigs$contig)]
ss.length = sapply(split(cb$length,cb$bin), sum)
ss.count = sapply(split(cb$length,cb$bin), length)
result = data.frame(bin=bins,
                    count=ss.count[match(bins,names(ss.count))],
                    length=ss.length[match(bins,names(ss.length))])

write.table(result, ofn, quote=F, row.names=F, sep="\t")
