library(stringr)
library(fst)


predict <- function(object, ...) {
    UseMethod("predict")
}

# This function generates predicted values for the output of BS.bigdata.
#' @param object An object of type BS.bigdata generated by smooth.bigdata.
#' @param new_data_set A big_data_set object generated by the big_data_set() function in big_data_set.R. This one will contain the list of Ft's, and is not expected to contain any Yt values.
#' @param old_data_set A big_data_set object generated by the big_data_set() function in big_data_set.R. The one used for the FFBS, specifically the FF procedure. Used to assist in generating values for Vt and Jt.
#' @param Vt left-covariance matrix of the noise parameter of Y.
#' @param Jt The covariance function for the new dataset being predicted and the old dataset used to generate the parameters.
#' @param out.head The file path to output the Y predictions. Defaults to the current directory with "Ypred_" as the header.
#' @param out.format The format of the files to read and write. Supports 'csv' and 'fst'.
#' @param fst.head The header of the name for the columns of both Y_t and F_t. Required to parse the columns from these files. Defaults to 'V' (i.e. the columns of both Y_t and F_t are labeled V1, V2, V3,...)
#' @param times The times to run the prediction function for. Takes in a select list of integers over which to generate predictions. Defaults to all times, i.e. 1 to the max time (length(object$st_files)).
#' @returns A list of files denoting the values of the predicted Y

predict.BS.bigdata <- function(object, new_data_set, old_data_set,
			       Vt_new = function(chunk) diag(nrow(chunk)),
			       Vt_old = function(chunk) diag(nrow(chunk)),
			       Jt = function(new_chunk, old_chunk) 
				       diag(max(nrow(new_chunk), nrow(old_chunk)))[1:nrow(old_chunk),1:nrow(new_chunk)], 
			       out.head = "./Ypred_", #grid.traversal.mode = "rowsnake",
			       out.format = "csv",
			       fst.head = "V",
			       times = seq(1, length(object$st_files)),
			       verbose = FALSE) {

    T <- length(object$st_files)
    # read in the columns of S_t. It has p(p+1)/2 columns.
#    p <- as.integer(system(sprintf("awk -F, '{print NF; exit}' %s", 
#					  object$St_files[1]), intern=TRUE))
    # solve quadratic formula
#    ptilde <- as.integer((-1 + sqrt(1 + 8 * p))/2)

    F_file_old <- old_data_set$get("F1")
    ptilde <- F_file_old$ncolblock

    F_file_new <- new_data_set$get("F1")
    file_name_F <- F_file_new$file_name
    p <- F_file_new$ncol
#    p <- as.integer(system(sprintf("awk -F, '{print NF; exit}' %s", file_name_F), intern=TRUE))
    
    Yt_pred_fname_list <- c()
    var_Yt_pred_fname_list <- c()

    for (t in times) {
    #for (t in 1:T) {
        if (out.format == "csv") {
            st_file <- read.csv(object$st_files[t])
            St_file <- read.csv(object$St_files[t])
	} else if (out.format == "fst") {
            st_file <- read.fst(object$st_files[t])
            St_file <- read.fst(object$St_files[t])
	}
	K <- nrow(st_file) - 1
       
	Y_file_old <- old_data_set$get(paste0("Y", t))
	file_name_Y <- Y_file_old$file_name

        F_file_new <- new_data_set$get(paste0("F", t))
        file_name_F_new <- F_file_new$file_name

	fnrow_new <- F_file_new$nrow
#        fnrow_new <- as.integer(str_split_1(system(sprintf("wc -l %s", file_name_F_new), intern=TRUE), ' ')[1]) -
#                                as.integer(F_file_new$header)
        # read in r and c from new_data_set
        r_new <- min(fnrow_new, F_file_new$nrowblock)
	
        c_F <- min(p, F_file_new$ncolblock)
        
        F_file_old <- old_data_set$get(paste0("F", t))
        file_name_F_old <- F_file_old$file_name

	fnrow_old <- F_file_old$nrow
#        fnrow_old <- as.integer(str_split_1(system(sprintf("wc -l %s", file_name_F_old), intern=TRUE), ' ')[1]) -
#                                as.integer(F_file_old$header)
	r_old <- min(fnrow_old, F_file_old$nrowblock)

	fncol <- Y_file_old$ncol
#        fncol <- as.integer(system(sprintf("awk -F, '{print NF; exit}' %s", file_name_Y), intern=TRUE))
        c <- min(fncol, Y_file_old$ncolblock)

        Y_grid <- Y_file_old$grid
#		generate_grid(fnrow_old, fncol, r_old, c, traversal.mode = grid.traversal.mode)

        # determine if column (row?) splitting is needed for our data. TODO: call it cell-splitting instead.
        split_col_F <- new_data_set$get("split_col_F")
        # only if split_col_F
        F_grid_new <- NULL
        if (split_col_F) {
            F_grid_new <- F_file_new$grid #generate_grid(fnrow_new, p, r_new, c_F, traversal.mode = grid.traversal.mode)
        }
    
        split_col_F <- old_data_set$get("split_col_F")
        # only if split_col_F
        F_grid_old <- NULL
        if (split_col_F) {
            F_grid_old <- F_file_old$grid #generate_grid(fnrow_old, p, r_old, c_F, traversal.mode = grid.traversal.mode)
        }

	#TODO: better to rearrange it in one big grid instead? i.e. fnrow_new x fncol? 
	#	And then populate it accordingly?
        Yt_pred_list <- matrix(nrow=K, ncol=r_new * c, data=0)
        var_Yt_pred_list <- matrix(nrow=K, ncol=r_new * (r_new + 1)/2, data=0)

	Yt_pred_fname <- paste0(out.head, sprintf("%d.%s", t, out.format))
        Yt_pred_fname_list <- c(Yt_pred_fname_list, Yt_pred_fname)
        
	var_Yt_pred_fname <- paste0(out.head, sprintf("_var_%d.%s", t, out.format))
	var_Yt_pred_fname_list <- c(var_Yt_pred_fname_list, var_Yt_pred_fname)

	for (k in 1:K) {
	    if (verbose) cat(sprintf("Predicting time %d partition %d/%d...", t, k, K))
	    # read in Ft from new_data_set
            Ft_new <- NULL
            if (is.null(F_grid_new)) {
                Ft_new <- if (out.format == "csv") {
			      as.matrix(read_big_csv(file_name_F_new, rows = ud, has_header=F_file_new$header,
					 ignore_header=TRUE))
                          } else if (out.format == "fst") {
                              as.matrix(read_big_fst(file_name_F_new, rows = ud)) 
                          }
            } else {
                ud_F <- unlist(F_grid_new[F_grid_new$block == k, c("U", "D")])
	        lr_F <- unlist(F_grid_new[F_grid_new$block == k, c("L", "R")])

                Ft_new <- if (out.format == "csv") {
                              as.matrix(read_big_csv(file_name_F_new, rows = ud_F,
                                        cols = unlist(F_grid_new[F_grid_new$block == k, c("L", "R")]),
                                        has_header=F_file_new$header,
			                ignore_header=TRUE))
		          } else if (out.format == "fst") {
                              as.matrix(read_big_fst(file_name_F_new, rows = ud_F,
						     cols = paste0(fst.head,
						            seq(lr_F[1], lr_F[2]))))
		          }
            }

            Ft_old <- NULL
            if (is.null(F_grid_old)) {
                Ft_old <- if (out.format == "csv") {
			      as.matrix(read_big_csv(file_name_F_old, rows = ud, has_header=F_file_old$header,
					 ignore_header=TRUE))
                          } else if (out.format == "fst") {
                              as.matrix(read_big_fst(file_name_F_old, rows = ud)) 
                          }
            } else {
                ud_F <- unlist(F_grid_old[F_grid_old$block == k, c("U", "D")])
	        lr_F <- unlist(F_grid_old[F_grid_old$block == k, c("L", "R")])
                Ft_old <- if (out.format == "csv") {
                              as.matrix(read_big_csv(file_name_F_old, rows = ud_F,
                                        cols = unlist(F_grid_old[F_grid_old$block == k, c("L", "R")]),
                                        has_header=F_file_old$header,
			                ignore_header=TRUE))
		          } else if (out.format == "fst") {
                              as.matrix(read_big_fst(file_name_F_old, rows = ud_F,
						     cols = paste0(fst.head,
							    seq(lr_F[1], lr_F[2]))))
		          }
            }
            
            s_t <- matrix(unlist(st_file[k,]), nrow=ptilde, ncol=c)
            S_t <- tri.to.sym(unlist(St_file[k,]), ptilde, lower=FALSE)
	    uS_t <- chol(S_t)


            V_t_new <- Vt_new(Ft_new)
            V_t_old <- Vt_old(Ft_old)
            uVt_old <- chol(V_t_old)

	    J_t <- Jt(Ft_new, Ft_old)

            ud <- unlist(Y_grid[Y_grid$block == k, c("U", "D")])
            lr <- unlist(Y_grid[Y_grid$block == k, c("L", "R")])

            # read in the Y
            Y_t <- if (out.format == "csv") {
		       as.matrix(read_big_csv(file_name_Y, rows = ud,
                                 cols = unlist(Y_grid[Y_grid$block == k, c("L", "R")]),
                                 has_header=Y_file_old$header,
			         ignore_header=TRUE))
	           } else if (out.format == "fst") {
		       as.matrix(read_big_fst(file_name_Y, rows = ud,
                                 cols = paste0(fst.head, seq(lr[1], lr[2]))))
	           }

            # build the predicted Y's. Also generate Jt and Vt here.
	    Yt_pred <- Ft_new %*% s_t + crossprod(J_t, backsolve(uVt_old, backsolve(uVt_old, Y_t - Ft_old %*% s_t, transpose=TRUE)))

            Ypost_var_old <- tcrossprod(Ft_old, uS_t)
            Ypost_var_old <- tcrossprod(Ypost_var_old) + V_t_old

	    uYpost_var_old <- chol(Ypost_var_old)
    
            Ypost_var_new <- tcrossprod(Ft_new, uS_t)
            Ypost_var_new <- tcrossprod(Ypost_var_new) + V_t_new
            Ypost_covar <- tcrossprod(Ft_old, uS_t) %*% tcrossprod(uS_t, Ft_new) + J_t
  
            uvaroldinv_covar <- backsolve(uYpost_var_old, Ypost_covar, transpose=TRUE)

            rowcov_Yt_pred <- Ypost_var_new - crossprod(uvaroldinv_covar)

            Yt_pred_list[k,] <- c(Yt_pred)
            var_Yt_pred_list[k,] <- rowcov_Yt_pred[upper.tri(rowcov_Yt_pred, diag=TRUE)]
            if (verbose) cat(sprintf("Prediction at time %d and partition %d/%d completed.\n", t, k, K))
	}
	# output recorded matrices to file.
	if (out.format == "csv") {
            write.csv(Yt_pred_list, file=Yt_pred_fname, row.names=FALSE)
            write.csv(var_Yt_pred_list, file=var_Yt_pred_fname, row.names=FALSE)
	} else if (out.format == "fst") {
            write.fst(as.data.frame(Yt_pred_list), Yt_pred_fname)
            write.fst(as.data.frame(var_Yt_pred_list), var_Yt_pred_fname)
	}
    }
    return(list(Yt_pred_means = Yt_pred_fname_list, 
		Yt_pred_rowcovs = var_Yt_pred_fname_list))
}

