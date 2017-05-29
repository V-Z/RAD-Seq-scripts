# R script to read standard input (single column of numbers) and calculate statistics
doc <- scan(file="stdin", quiet=TRUE)
summary(doc)
# cat(min(doc), max(doc), median(doc), mean(doc), sep="\t")
