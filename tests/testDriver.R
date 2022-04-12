library("testit")
icon = list(
	success = "\U2713",
	fail = "\U2717"
)

message("Testing Basic intervals tasks")
source("src/intervals.R")

message("Testing countOverlapsSimple")

set.seed(1)
starts = floor(runif(1000000, 100, 5000))
lengths = abs(floor(rnorm(1000000, 100, 25))+10)
startsOrdered = starts[order(starts)]
endsOrdered = startsOrdered + lengths
q = data.frame(start=2200, end=2250)
df = data.frame(starts=startsOrdered, ends=endsOrdered)

ol = countOverlapsSimple(q, df)

assert("countOverlapsSimple", {
	(ol == 32789)
})

message(icon$success, " Passed countOverlapsSimple basic tests!")

suppressMessages(suppressWarnings(library("GenomicRanges")))

# Now we'll generate some random data to ensure that your 
# countOverlapsSimple function is working
for (i in 1:50) {
	randdb = randomGRanges()
	fixedquery = list(start=1000000, end=2000000)
	qgr = GRanges("chrX", ranges=IRanges(fixedquery$start, fixedquery$end))
	dbdf = data.frame(starts=start(randdb), ends=end(randdb))
	cos = countOverlapsSimple(fixedquery, dbdf)
	check = GenomicRanges::countOverlaps(qgr, randdb)
	assert("countOverlapsSimple Random Tests", {
		(cos == check)
	})
}

message(icon$success, " Passed countOverlapsSimple random tests!")
message("Testing calculateJaccardScore")

filesList = c(
	"data/helas3_ctcf.narrowPeak.gz",
	"data/helas3_jun.narrowPeak.gz",
	"data/hepg2_ctcf.narrowPeak.gz",
	"data/hepg2_jun.narrowPeak.gz"
)

files = lapply(filesList, read.table)
granges = lapply(files, function(x) {
	GRanges(seqnames=x$V1, ranges=IRanges(x$V2, x$V3))
})

assert("calculateJaccardScore", {
	(abs(calculateJaccardScore(granges[[1]], granges[[2]]) - 0.02053016) < .005)
	(abs(calculateJaccardScore(granges[[1]], granges[[3]]) - 0.607247) < .005)
	(abs(calculateJaccardScore(granges[[1]], granges[[4]]) - 0.023433) < .005)
	(abs(calculateJaccardScore(granges[[2]], granges[[3]]) - 0.01328) < .005)
	(abs(calculateJaccardScore(granges[[2]], granges[[4]]) - 0.1665181) < .005)
})

message(icon$success, " Passed calculateJaccardScore tests!")

message("Testing pairwiseJaccard")

pwj = pairwiseJaccard(granges)

message(icon$success, " Passed pairwiseJaccard tests!")

message(icon$success, " All Interval Basics tests passed!")
