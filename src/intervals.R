# Generates a random GRanges object
#
# You can use this function to get a random interval list for testing.
# Don't change this function, it's just for you to use while
# you work on the other functions
#' @examples
# gr1 = randomGRanges()
randomGRanges = function(){
	starts = sample(1:1e7, 500)
	sizes = sample(1:1e4, 500)
	ends = starts + sizes
	ret = GenomicRanges::GRanges(seqnames="chrX", ranges=IRanges::IRanges(starts,ends))
	GenomicRanges::reduce(ret)
}

# Counts overlaps between query and database
#
# This function counts the number of intervals from a 
# database set of intervals that overlap with a single
# query interval. 
#
#' @param query  
#'   Interval to check for overlaps. Should be a list 
#'   object with 'start' and 'end' values.
#' @param database
#'   A data.frame, with rows corresponding to intervals,
#'   with 'starts' and 'ends' as columns.
#' @return
#'   Number of overlaps counted (as a numeric)
#' @examples
#'   query = list(start=50, end=60)
#'   database = data.frame(starts=c(35, 45, 55), ends=c(45, 65, 85))
#'   countOverlapsSimple(query, database)  # returns 2
countOverlapsSimple = function(query, database) {
  
  # wrapper for function overlapping
  .doesOverlap = function(A, B) {
    # map to better names
    A_start = A[1]
    A_end = A[2]
    B_start = B$start
    B_end = B$end
    
    # check olap
    if (A_start <= B_end && A_end >= B_start) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  }
  
  # calculate all overlaps
  database$overlaps_query <- apply(database, 1, .doesOverlap, B=query)
  
  # debug
  # print(sum(database$overlaps_query, na.rm=TRUE))
  
  # return count
  return(sum(database$overlaps_query, na.rm=TRUE))

}

# Measure the Jaccard similarity between two interval sets
#
# This function calculates the Jaccard Score between two
# given interval sets, provided as GenomicRanges::GRanges 
# objects. The Jaccard score is computed as the intersection
# between the two sets divided by the union of the two sets.
# @param gr1  First GRanges object
# @param gr2  Second GRanges object
# @return The Jaccard score (as numeric)
calculateJaccardScore = function(gr1, gr2){
	# Implement this function
  return(
    as.numeric(
      length(intersect(gr1, gr2))/length(union(gr1, gr2))
    )
  )
}

# Calculate pairwise Jaccard similarity among several interval sets
#
# This function makes use of \code{calculateJaccardScore}. It simply
# loops through each pairwise comparison and calculates.
#' Round the result to 3 significant figures using \code{signif}
#' @param lst A base R list of GRanges objects to compare
#' @return
#'   A matrix of size n-by-n, where n is the number of elements in lst.
#' @examples
#' lst = replicate(10, randomGRanges())
#' pairwiseJaccard(lst)
pairwiseJaccard = function(lst) {
	# Compute the pairwise Jaccard Score and return in matrix form
  N = length(lst)
  J_MAT = matrix(nrow=N, ncol=N)
  
  # loop through
  for (i in 1:N) {
    for (j in 1:N) {
      # Return only 3 significant figures
      J_MAT[i,j] = signif(calculateJaccardScore(lst[[i]], lst[[j]]), digits=3)
    }
  }
  
  return(J_MAT)
}



