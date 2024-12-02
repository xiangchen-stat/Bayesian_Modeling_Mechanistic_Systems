library(MNIW)
fnrow <- 1000
fncol <- 1000
bnrow <- 50
bncol <- 50

mygrid <- generate_grid(fnrow, fncol, bnrow, bncol, 
                        traversal.mode = "rowsnake", is.flexible = FALSE)

# Test that cell 1 is equal to c(1,500,1,500)
iv <- unlist(mygrid[mygrid$block_no == 1, c("U","D","L","R")])
names(iv) <- NULL
expect_equal(iv, c(1,500,1,500))

# Test that cell 20 is equal to c(1,500, 9501, 10000)
iv <- unlist(mygrid[mygrid$block_no == 20, c("U","D","L","R")])
names(iv) <- NULL
expect_equal(iv, c(1,500,9501,10000))

# Test that cell 21 is equal to c(501, 1000, 9501, 10000)
iv <- unlist(mygrid[mygrid$block_no == 21, c("U","D","L","R")])
names(iv) <- NULL
expect_equal(iv, c(501,1000,9501,10000))

# Test that cell 22 is equal to c(501, 1000, 9001, 9500)
iv <- unlist(mygrid[mygrid$block_no == 22, c("U","D","L","R")])
names(iv) <- NULL
expect_equal(iv, c(501,1000,9001,9500))

# Test that cell 400 is equal to c(9501, 10000, 9501, 10000)
iv <- unlist(mygrid[mygrid$block_no == 400, c("U","D","L","R")])
names(iv) <- NULL
expect_equal(iv, c(9501,10000,1,500))