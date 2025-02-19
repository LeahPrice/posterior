#' @importFrom rlang eval_tidy is_missing missing_arg dots_list
#' @export
`[[.rvar` <- function(x, i, ...) {
  index <- check_rvar_yank_index(x, i, ...)

  if (length(index) == 1) {
    # single element selection => collapse the dims so we can select directly using i
    .dim = dim(x)
    if (length(.dim) != 1) {
      # we only collapse dims if necessary since this will drop dimnames (which
      # would prevent single-element by-name selection for 1d rvars)
      dim(x) <- prod(.dim)
    }
    .draws <- draws_of(x)[, i, drop = FALSE]
    dimnames(.draws) <- NULL
    out <- new_rvar(.draws, .nchains = nchains(x))
  } else if (length(index) == length(dim(x))) {
    # multiple element selection => must have exactly the right number of dims
    .draws <- eval_tidy(expr(draws_of(x)[, !!!index, drop = FALSE]))
    # must do drop manually in case the draws dimension has only 1 draw
    dim(.draws) <- c(ndraws(x), 1)
    out <- new_rvar(.draws, .nchains = nchains(x))
  } else {
    stop_no_call("subscript out of bounds")
  }
  out
}

#' @export
`[[<-.rvar` <- function(x, i, ..., value) {
  value <- vec_cast(value, x)
  c(x, value) %<-% conform_rvar_ndraws_nchains(list(x, value))
  value <- check_rvar_dims_first(value, new_rvar(0))
  index <- check_rvar_yank_index(x, i, ...)

  if (length(index) == 1) {
    .dim = dim(x)

    if (length(.dim) == 1 && i > length(x)) {
      # unidimensional indexing allows array extension; extend the array
      # then do the assignment
      x <- x[seq_len(max(i, na.rm = TRUE))]
      draws_of(x)[, i] <- draws_of(value)
    } else {
      # single element selection => collapse the dims so we can select directly using i
      .dimnames = dimnames(draws_of(x)) # to restore later
      if (length(.dim) != 1) {
        # we only collapse dims if necessary since this will drop dimnames (which
        # would prevent single-element by-name selection for 1d rvars)
        dim(x) <- prod(.dim)
      }
      draws_of(x)[, i] <- draws_of(value)
      dim(x) <- .dim
      dimnames(draws_of(x)) <- .dimnames
    }
  } else if (length(index) == length(dim(x))) {
    # multiple element selection => must have exactly the right number of dims
    x <- eval_tidy(expr({
      draws_of(x)[, !!!index] <- draws_of(value)
      x
    }))
  } else {
    stop_no_call("subscript out of bounds")
  }
  x
}

#' @importFrom rlang dots_list eval_tidy is_missing missing_arg expr
#' @export
`[.rvar` <- function(x, ..., drop = FALSE) {
  check_rvar_subset_indices(x, ...)
  .draws = draws_of(x)
  .dim = dim(.draws)

  # clean up the indices: because we have to index using the multi-index
  # notation x[,...] (to account for draws) and this notation has slightly
  # different semantics for NAs from the way that x[i] works, have to do a bit
  # of translation here
  index = dots_list(..., .ignore_empty = "none", .preserve_empty = TRUE)
  for (i in seq_along(index)) {
    if (is.numeric(index[[i]])) {
      # numeric indices outside the range of the corresponding dimension
      # should create NAs; but array indexing doesn't do this (it throws
      # an error), so we adjust the indices to do so.

      dim_i_length <- if (i == 1 && length(index) == 1) {
        # for x[i] style indexing of multidimensional arrays we will flatten
        # the array before indexing, so the max of the dim length will be
        # the length of x
        length(x)
      } else {
        .dim[[i + 1]]
      }

      index[[i]][index[[i]] > dim_i_length] <- NA_integer_
    }
  }

  if (length(index) == 1) {
    # indexing by a single dimension, (could be numerical indexing along one
    # dimension --- even for multidimensional arrays; logical indexing; or matrix indexing)

    if (is.matrix(index[[1]]) && ncol(index[[1]]) == length(.dim) - 1) {
      # matrix-based indexing, like x[cbind(2,1,2)]
      # => translate matrix-based indices into unidimensional indices
      index[[1]] <- matrix_to_index(index[[1]], .dim[-1])
    }

    if (is_missing(index[[1]])) {
      # if we only have one index and it is missing entirely, the call must
      # have been for x[], which actually means no index at all
      index <- list()
    } else if (length(dim(.draws)) > 2) {
      # draws have > 2 dims => array has > 1 dimension => must flatten array so
      # that the index can be applied along all the elements of the array
      dim(.draws) <- c(dim(.draws)[[1]], length(x))
    }
  } else if (length(index) < length(dim(.draws)) - 1) {
    # fill in final indices with missing arguments
    index[seq(length(index) + 1, length(dim(.draws)) - 1)] = list(missing_arg())
  }

  x <- eval_tidy(expr(
    new_rvar(.draws[, !!!index, drop = FALSE], .nchains = nchains(x))
  ))

  if (drop) {
    x <- drop(x)
  }

  x
}

#' @export
`[<-.rvar` <- function(x, i, ..., value) {
  if (missing(i)) i = missing_arg()
  if (length(dim(x)) == 1 && !missing(i) && any(i > length(x), na.rm = TRUE)) {
    # unidimensional indexing allows array extension; extend the array
    # before we do the assignment
    x <- x[seq_len(max(i, na.rm = TRUE))]
  }

  value <- vec_cast(value, x)
  c(x, value) %<-% conform_rvar_ndraws_nchains(list(x, value))

  if (missing(...)) {
    # index over entire array: flatten array so that the  index
    # can be applied along all the elements of the array, then invert after assignment
    original_dim <- dim(draws_of(x))
    original_dimnames <- dimnames(draws_of(x))

    if (is.matrix(i) && ncol(i) == length(original_dim) - 1) {
      # matrix-based indexing, like x[cbind(2,1,2)] <- y
      # => translate matrix-based indices into unidimensional indices
      i <- matrix_to_index(i, original_dim[-1])
    }

    #flatten and assign
    dim(x) <- length(x)
    draws_of(x)[,i] <- draws_of(value)

    # unflatten and restore dimnames
    dim(draws_of(x)) <- original_dim
    dimnames(draws_of(x)) <- original_dimnames
  } else {
    draws_of(x)[,i,...] <- draws_of(value)
  }

  x
}



# slicing helpers ---------------------------------------------------------

# Given m, a matrix with ncol(m) == length(dim), where each row specifies the
# index of a single cell from an array of shape == dim, return a vector of
# length nrow(m) giving the unidimensional indices of that array corresponding
# to the cells specified by each row in m
# e.g. if we have an array x with dim(x) = c(2,3,4), we might want to
# translate an index of style x[cbind(2,1,2)] into x[i] so that we can index
# index x using the single index i instead of cbind(2,1,2).
# matrix_to_index(cbind(2,1,2), dim(x)) does that translation.
matrix_to_index <- function(m, dim) {
  cumdim <- cumprod(c(1, dim[-length(dim)]))
  as.vector((m - 1) %*% cumdim + 1)
}
