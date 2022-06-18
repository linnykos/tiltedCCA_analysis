input_obj <- eSVD_obj
latest_Fit <- eSVD2:::.get_object(eSVD_obj = input_obj,
                          what_obj = "latest_Fit",
                          which_fit = NULL)
covariates <- eSVD2:::.get_object(eSVD_obj = input_obj,
                          what_obj = "covariates",
                          which_fit = NULL)
x_mat <- eSVD2:::.get_object(eSVD_obj = input_obj,
                     what_obj = "x_mat",
                     which_fit = latest_Fit)
y_mat <- eSVD2:::.get_object(eSVD_obj = input_obj,
                     what_obj = "y_mat",
                     which_fit = latest_Fit)
z_mat <- eSVD2:::.get_object(eSVD_obj = input_obj,
                     what_obj = "z_mat",
                     which_fit = latest_Fit)
nat_mat1 <- tcrossprod(x_mat, y_mat)
nat_mat2 <- tcrossprod(covariates, z_mat)

variable = "RPS24"
vec <- exp(nat_mat1[,variable] + nat_mat2[,variable])
quantile(vec)

quantile(as.numeric(input_obj$dat[,variable]))
