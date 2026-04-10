options(stringsAsFactors=FALSE)
exprDir <- "./"

samples <- sub("_filt_v_genome_sort_rs_count.txt$", "", 
               list.files(exprDir, pattern="_filt_v_genome_sort_rs_count.txt$"))
qexprList <- lapply(samples, function(s) {
    read.delim(sprintf("%s/%s_filt_v_genome_sort_rs_count.txt", exprDir, s), 
               header=FALSE, row.names=1, na.strings=c("N/A"))
})
names(qexprList) <- samples

qeTable <- do.call(cbind.data.frame, lapply(qexprList, function(qe) { qe[[2]] }))
colnames(qeTable) <- samples
row.names(qeTable) <- row.names(qexprList[[1]])

write.csv(qeTable, sprintf("%s/exprtable-all.csv", exprDir))

qeSummary <- data.frame(riboMean = rowMeans(qeTable),
                        row.names=row.names(qeTable))
qeSummary$riboBucket = 2**floor(log2(qeSummary$riboMean))

write.csv(qeTable, sprintf("%s/exprtable-summary.csv", exprDir))
write.table(table(qeSummary$riboBucket), 
            sprintf("%s/exprtable-binned-summary.txt", exprDir), sep="\t")

present <- qeTable[qeSummary$riboBucket > 31,]
write.csv(present, sprintf("%s/exprtable-present.csv", exprDir))
cor(qeTable[qeSummary$riboBucket > 31,], method="spearman")
