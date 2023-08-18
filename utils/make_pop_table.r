#!/bin/env Rscript

options(warn=1)

# get script name
all.args = commandArgs(F)
fn.arg = "--file="
script.name = sub(fn.arg, "", all.args[grep(fn.arg, all.args)])

args = commandArgs(T)
if (length(args) == 0) {
  cat(sprintf("usage: %s <base dir> <filename> <ofn> <sample1> <sample2> ...\n", script.name))
  q(status=1)
}

dir = args[1]
sfn = args[2]
ofn = args[3]
ids = args[4:length(args)]

fns = paste0(dir, "/", ids, "/", sfn)
result = data.frame(id=ids, fn=fns)

write.table(result, ofn, quote=F, row.names=F, sep="\t")
