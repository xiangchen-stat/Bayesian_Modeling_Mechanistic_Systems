make_seed_obj <- function(seed_list = NULL, seed_ix = 1,
			seed_now = 0, seed_step = 1, seedarr_explicit = FALSE) {
    if (is.null(seed_now)) {
        out <- list("seed_now" = NULL)
        return(out)
    }

    out <- if (seedarr_explicit) {
               list("seeds"=seed_list,"seed_ix"=seed_ix, "seed_step"=seed_step,
		    "seedarr_explicit"=seedarr_explicit)
           } else {
               list("seed_now"=seed_now, "seed_step"=seed_step,
		    "seedarr_explicit"=seedarr_explicit)
           }
    class(out) <- "seed_obj"
    return(out)
    
}

#' Updates and sets the seed according to preset values.
#' @param seed_obj A seed object created by make_seed_obj()
#' @param value The number to increment the seed index by if the list of seeds is given explicitly, or the current seed if only the current seed and the amount is given.

`update_seed<-` <- function(seed_obj, value = 1) {
    if (is.null(seed_obj)) {
        set.seed(NULL)
        return(NULL)
    }
    if (seed_obj$seedarr_explicit) {
        seed_obj$seed_ix <- seed_obj$seed_ix + value
        set.seed(seed_obj$seeds[seed_obj$seed_ix])
    } else {
        seed_obj$seed_now <- seed_obj$seed_now + value
        set.seed(seed_obj$seed_now)
    }
    seed_obj
}	

#' Updates and sets the seed according to preset values.
#' @param seed_obj A seed object created by make_seed_obj()
#' @param value The number to increment the seed index by if the list of seeds is given explicitly, or the current seed if only the current seed and the amount is given.

update_seed_newcopy <- function(seed_obj, value = 1) {
    if (is.null(seed_obj)) {
        set.seed(NULL)
        return(NULL)
    }
    if (seed_obj$seedarr_explicit) {
        seed_obj$seed_ix <- seed_obj$seed_ix + value
        set.seed(seed_obj$seeds[seed_obj$seed_ix])
    } else {
        seed_obj$seed_now <- seed_obj$seed_now + value
        set.seed(seed_obj$seed_now)
    }
    return(seed_obj)
}	


