library(readr)
library(fst)
library(collections)
library(stringr)

options(readr.show_progress = FALSE)

#' This function takes in a path to a specific file and the number of rows and columns to read in at one time.
#' @param filename The name of the file to read in
#' @param fnrow The number of rows in your file.
#' @param fncol The number of columns in your file.
#' @param nrowblock The number of rows to read in one block. Defaults to infinity, in which all rows are read.
#' @param ncolblock The number of columns to read in one block. Defaults to infinity, in which all columns are read.
#' @param header Whether or not to read the header of the file. Defaults to TRUE.
#' @returns A big_data_file object containing the input parameters.

big_data_file <- function(filename, fnrow, fncol, nrowblock = Inf, ncolblock = Inf, header=TRUE, grid.traversal.mode=NULL) {
    #TODO: read first line of filename to store headers into the object?
    # speed bottleneck here. Need a faster way to compute this.	
    # Record number of rows and columns the file has
#    fnrow <- as.integer(str_split_1(system(sprintf("wc -l %s", filename), intern=TRUE), ' ')[1]) -
#                 as.integer(header)
#    fncol <- as.integer(system(sprintf("awk -F, '{print NF; exit}' %s", filename), intern=TRUE))

    bnrow <- min(fnrow, nrowblock)
    bncol <- min(fncol, ncolblock)

    file_grid <- NULL
    K <- NULL
    if (!is.null(grid.traversal.mode)) {
        file_grid <- generate_grid(fnrow, fncol, bnrow, bncol, 
				   traversal.mode=grid.traversal.mode)
        K <- nrow(file_grid)
    }

    out <- list(file_name = filename,
	        nrow = fnrow, ncol = fncol,	
		nrowblock = nrowblock, ncolblock = ncolblock, 
		header=header,
                grid=file_grid,
                K=K)
    class(out) <- "big_data_file"
    return(out)
}

#' Reads part of a big fst file in. Analogous to read_big_csv, but for fst files.
#' @param filename The file to read in.
#' @param rows An ordered pair denoting the first and last row to read. Defaults to c(1,Inf), in which the function reads in all rows.
#' @param cols A list of character names for the columns to read in. If NULL, will read in all columns. Default: NULL.
#' @param ignore_header Whether or not to ignore the header when outputting the file. Defaults to FALSE.

read_big_fst <- function(filename, rows=c(1,Inf), cols=NULL, ignore_header=FALSE) {
    from = rows[1]
    to = rows[2]
    if (!is.finite(to)) to <- NULL

    out <- read_fst(filename, columns=cols, from=from, to=to)
    return(out)
}

#' This function reads in part of a big data file and outputs the CSV as a data.frame.
#' @param filename The file to read in.
#' @param rows The rows to read in. An ordered pair denoting the first and last row to read. Defaults to c(1,Inf), in which the function reads in all rows.
#' @param cols The columns to read in. An ordered pair denoting the first and last columns to read. Defaults to c(1,Inf), in which the function reads in all columns.
#' @param has_header Whether or not the file has a header to be read in. Defaults to TRUE.
#' @param ignore_header Whether or not to ignore the header when outputting the file. Defaults to FALSE.

read_big_csv <- function(filename, rows = c(1,Inf), cols = c(1,Inf),
			 has_header=TRUE, ignore_header=FALSE) {

    if (is.finite(cols[2])) {
        if (!ignore_header) headers <- as.character(as.vector(read_csv(filename,
					       col_names = FALSE, n_max = 1,
					       col_select = cols[1]:cols[2],
					       show_col_types=FALSE)))
        out <- read_csv(file=filename, col_names = FALSE,
                   skip = rows[1] - 1 + as.integer(has_header), 
		   n_max = rows[2] - rows[1] + 1,
                   col_select = cols[1]:cols[2],
                   show_col_types=FALSE)

    } else {
        #TODO: something is going on with this read_csv that is throwing an error, but when I run it in the R shell, it works fine. Oddly, neglecting this option does not throw an error. 
        if (!ignore_header) headers <- read_csv(filename, col_names = FALSE, n_max = 1,
               show_col_types=FALSE)
        out <- read_csv(file=filename, col_names = FALSE,
                   skip = rows[1] - 1 + as.integer(has_header), 
		   n_max = rows[2] - rows[1] + 1,
                   show_col_types=FALSE)
    }

    if (has_header & !ignore_header) colnames(out) <- headers
    return(out)
}


#' This function encapsulates a set of file names into big_data_file objects and then groups them into a numerically-indexed list for the FFBS to access its elements more easily.
#' @param Y_filename_list The list of file names to read in.
#' @param F_filename_list The file containing the covariates F, one per Y.
#' @param fnrow_list A list of the number of rows in each outcome file Y.
#' @param fncol_list A list of the number of columns in each outcome file Y.
#' @param F_fncol_list A list of the number of columns in each covariate file F.
#' @param nrowblock_list The number of rows to read in one block. The big_data_set expects one positive integer per filename in Y_filename_list. Defaults to infinity for all.
#' @param ncolblock The number of columns to read in one block. The big_data_set expects one positive integer per filename in Y_filename_list. Defaults to infinity for all.
#' @param F_ncolblock_list The number of columns to read in one block. Defaults to ncolblock_list.
#' @param header Whether or not to read the header of the file. Defaults to TRUE.
#' @param F_header Whether or not to read the header of the covariates F. Defaults to the header option.
#' @param split_col_F Whether or not to split the columns of F into blocks. Defaults to FALSE.
#' @param grid.traversal.mode The traversal procedure for the data blocks for both Y and F. Options include "rowsnake" and "lr". Must be specified.
#' @param uniform_grid Whether or not to use the same block size for all files for both Y and F. If TRUE, then it will take the first element of nrowblock_list and ncolblock_list and apply it to all Y's, same for F if split_col_F is set to TRUE. Defaults to TRUE.
#' @returns a big_data_set object containing a list of big_data_file objects. The object's underlying implementation is a dict(), so its elements are accessed using $get().

#TODO: add an option for whether or not the grid is exact??
big_data_set <- function(Y_filename_list,
			 F_filename_list,
			 fnrow_list,
			 fncol_list,
			 F_fncol_list,
			 nrowblock_list = rep(Inf, length(Y_filename_list)),
			 ncolblock_list = rep(Inf, length(Y_filename_list)),
			 F_ncolblock_list = ncolblock_list, 
			 header=TRUE,
			 F_header=header,
			 split_col_F = FALSE,
			 grid.traversal.mode = NULL,
			 uniform_grid = TRUE) {

    out = dict()
    out$set("split_col_F", split_col_F)
    out$set("uniform_grid", uniform_grid)

    nrowblock = nrowblock_list[1]
    ncolblock = ncolblock_list[1]

    F_ncolblock = F_ncolblock_list[1]

    T <- length(Y_filename_list)

    for (i in 1:T) {
        Y_ix <- sprintf("Y%d", i)
        F_ix <- sprintf("F%d", i)
        if (!uniform_grid) {
            nrowblock = nrowblock_list[i]
            ncolblock = ncolblock_list[i]
	    if (split_col_F) F_ncolblock = F_ncolblock_list[i]
	}

        out$set(Y_ix, big_data_file(Y_filename_list[i],
				 fnrow=fnrow_list[i],
				 fncol=fncol_list[i],
				 nrowblock=nrowblock,
				 ncolblock=ncolblock,
				 header=header,
				 grid.traversal.mode=grid.traversal.mode))
        # just store the file name. F's block size will be controlled by the block sizes of the Y's.
        out$set(F_ix, big_data_file(F_filename_list[i],
				 fnrow=fnrow_list[i],
				 fncol=F_fncol_list[i],
				 nrowblock=nrowblock,
				 ncolblock=F_ncolblock,
				 header=F_header,
				 grid.traversal.mode=grid.traversal.mode))
	#Assert that out$F(k)$K == out$Y(K)$K. Else stop and throw an error.
	if (split_col_F & (out$get(F_ix)$K != out$get(Y_ix)$K)) {
	    stop(sprintf("Unequal number of blocks for Y and F for file %d! Check your block lists.", i))
	}
    }

    sum_K <- 0
    for (i in 1:T) {
        Y_ix <- sprintf("Y%d", i)
        sum_K <- sum_K + out$get(Y_ix)$K
    }
    out$set("sum_K", sum_K)

    class(out) <- "big_data_set"
    return(out)
}

