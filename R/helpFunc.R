#' Return the first or last part of a list (from
#' https://gist.github.com/pimentel/256fc8c9b5191da63819)
#'
#' Returns the first or last part of a list. Instead of returning the first
#' n entries as the standard head() does, it attempts to call head()
#' recursively on the entries in the list. If it fails, it will return the
#' particular entry (standard behavior).
#' @param obj a list object
#' @param n a single integer. If positive, prints the first n items for the
#' list and all entries in the list. If negative, prints all but the last
#' n items in the list.
#' @return a list of length n, with items in the list of length n
#'
#' @export
headList = function(obj, n = 6L)
{
  stopifnot(length(n) == 1L)
  origN = n
  n = if (n < 0L)
    max(length(obj) + n, 0L)
  else min(n, length(obj))
  lapply(obj[seq_len(n)], function(x)
  {
    tryCatch({
      head(x, origN)
    }, error = function(e) {
      x
    })
  })
}
