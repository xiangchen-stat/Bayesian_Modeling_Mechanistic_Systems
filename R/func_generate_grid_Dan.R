#' Generate the indices for a grid for big data traversal. Indexes the grid by a number from 1 to K, where K is the total number of blocks partitioned using the grid.
#' @param fnrow The number of rows in the file being read.
#' @param fncol The number of columns in the file being read.
#' @param bnrow The number of rows in each block being read.
#' @param bncol The number of columns in each block being read.
#' @param traversal.mode The way to traverse the grid. "rowsnake" refers to alternating left-right traversal. "lr" refers to left-right traversal throughout. Defaults to "rowsnake".
#' @param is.flexible Whether or not to throw an error if bnrow or bncol do not perfectly divide fnrow and fncol respectively. If set, will override both is.flexible.row and is.flexible.col. Defaults to NULL.
#' @param is.flexible.row Whether or not to throw an error if bnrow does not perfectly divide fnrow. Defaults to FALSE.
#' @param is.flexible.col Whether or not to throw an error if bncol does not perfectly divide fncol. Defaults to FALSE.
#' @returns A data.frame with the columns for block number, row upper index, row lower index, column left index, and column right index, in that order.

#TODO: make this more flexible by allowing fnrow and fncol to be intervals of integers?
generate_grid <- function(fnrow, fncol, bnrow, bncol, traversal.mode = "rowsnake",
			  is.flexible=NULL, is.flexible.row = FALSE, is.flexible.col = FALSE) {
    if (!is.null(is.flexible)) {
        is.flexible.row <- is.flexible
        is.flexible.col <- is.flexible
    }
    if (!is.flexible.row & (fnrow %% bnrow != 0)) stop(sprintf("Number of rows per block %d must divide number of rows in file %d.", bnrow, fnrow))
    if (!is.flexible.col & (fncol %% bncol != 0)) stop(sprintf("Number of cols per block %d must divide number of cols in file %d.", bncol, fncol))

    if (!(is.flexible.row | is.flexible.col)) return(generate.grid.exact(fnrow, fncol, bnrow, bncol, traversal.mode))

    if (traversal.mode == "rowsnake") return(generate.grid.rowsnake(fnrow, fncol, bnrow, bncol))
    if (traversal.mode == "lr") return(generate.grid.lr(fnrow, fncol, bnrow, bncol))
}

# generate the grid analytically. Does not require a for loop and is designed to be fast.
generate.grid.exact <- function(fnrow, fncol, bnrow, bncol, traversal.mode) {
    Kr <- fnrow/bnrow
    Kc <- fncol/bncol
    K <- Kr * Kc

    out <- data.frame(block_no = seq(1,K), U = rep(seq(1, fnrow, bnrow), each=Kc), D = rep(seq(bnrow, fnrow, bnrow), each=Kc), 
		      L = rep(seq(1, fncol, bncol), Kr), R = rep(seq(bncol, fncol, bncol), Kr))

    if (traversal.mode == "rowsnake") {
        # flip the left-right of every even row (UD), i.e. where block_no is in a certain range.
        out[(out$D %% (2 * bnrow) == 0), "L"] <- rev(rep(seq(1, fncol, bncol), Kr/2))
        out[(out$D %% (2 * bnrow) == 0), "R"] <- rev(rep(seq(bncol, fncol, bncol), Kr/2))
    }

    return(out)
}

# The rowsnake version of the grid generation.
generate.grid.rowsnake <- function(fnrow, fncol, bnrow, bncol) {
    Kr <- ceiling(fnrow/bnrow)
    Kc <- ceiling(fncol/bncol)

    K <- Kr * Kc

    U <- 1
    D <- U + bnrow - 1
    L <- 1
    R <- L + bncol - 1

    out <- data.frame(block_no = c(1), U = c(U), D = c(D), L = c(L), R = c(R))
    # the direction to traverse
    #TODO: toggle in an option in the future on where to start.
    traversal <- "right"
    traversal.last <- "right"

    #TODO: There should be a way to analytically calculate each partition of the grid without 
    #		regenerating the dataframe each time or throwing this in a for loop.
    # 	Check by running testthat.
    for (k in 2:K) {
        if (traversal == "right") {
            L <- min(L + bncol, fncol)
            R <- min(R + bncol, fncol)
            # cannot go further
            if (R == fncol) {
	        traversal <- "down"
                traversal.last <- "right"        
	    }
        } else if (traversal == "down") {
            # update columns
            U <- min(U + bnrow, fnrow)
            D <- min(D + bnrow, fnrow)
            if (traversal.last == "right") traversal <- "left"
	    if (traversal.last == "left")  traversal <- "right"

        } else if (traversal == "left") {
            R <- ceiling((R - bncol)/bncol) * bncol
            L <- R - bncol + 1

            # cannot go further
            if (L == 1) {
	        traversal <- "down"
                traversal.last <- "left"        
	    }
        }
        out <- rbind(out, data.frame(block_no = c(k), U = c(U), D = c(D), L = c(L), R = c(R)))
    }
    return(out)
}

generate.grid.lr <- function(fnrow, fncol, bnrow, bncol) {
    Kr <- ceiling(fnrow/bnrow)
    Kc <- ceiling(fncol/bncol)

    K <- Kr * Kc

    U <- 1
    D <- U + bnrow - 1
    L <- 1
    R <- L + bncol - 1

    out <- data.frame(block_no = c(1), U = c(U), D = c(D), L = c(L), R = c(R))
    for (k in 2:K) {
	if (R == fncol) {
            U <- min(U + bnrow, fnrow)
            D <- min(D + bnrow, fnrow)
            L <- 1
	    R <- L + bncol - 1
	} else {
            L <- min(L + bncol, fncol)
            R <- min(R + bncol, fncol)
	}
        out <- rbind(out, data.frame(block_no = c(k), U = c(U), D = c(D), L = c(L), R = c(R)))
    }
    return(out)
}

