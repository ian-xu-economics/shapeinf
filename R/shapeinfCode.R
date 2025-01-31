#' Simulate a dataset from a specified data-generating process
#'
#' This function simulates data of size \code{n} from a user-defined data-generating process (DGP).
#' The DGP includes:
#' \itemize{
#'   \item A support and probability distribution for \code{Z}.
#'   \item A function to compute the propensity score for treatment assignment.
#'   \item Model parameters (mtrs) for generating the outcome \code{Y} conditional on treatment \code{D}.
#' }
#'
#' @param dgp A list describing the data-generating process. Must contain:
#'   \itemize{
#'     \item \code{suppZ}: Vector of possible values of \code{Z}.
#'     \item \code{densZ}: Probability distribution over \code{suppZ}.
#'     \item \code{find_pscore}: A function taking a value of \code{Z} and returning a propensity score.
#'     \item \code{mtrs}: A list of length 2, each entry containing:
#'       \itemize{
#'         \item \code{theta}: Numeric vector of parameters.
#'         \item \code{basis}: A list \code{b} of basis functions for the outcome model.
#'       }
#'   }
#' @param n An integer specifying the number of observations to simulate.
#' @param extra.details Logical; if TRUE, returns probability y equals 1 and u.
#' 
#' @return A \code{tibble} with columns:
#'   \itemize{
#'     \item \code{y}: Simulated outcome (0/1).
#'     \item \code{probability.y}: The probability of \code{y} being 1, given \code{D} and \code{Z}.
#'     \item \code{d}: Treatment assignment (0/1).
#'     \item \code{u}: The uniform random number used in drawing \code{D} and \code{Y}.
#'     \item \code{z}: The drawn value of the covariate.
#'   }
#'
#' @examples
#' # Suppose we have a simple DGP:
#' dgp_example <- list(
#'   suppZ = c(0, 1),
#'   densZ = c(0.5, 0.5),
#'   find_pscore = function(z) { ifelse(z == 0, 0.3, 0.7) },
#'   mtrs = list(
#'     list(
#'       theta = c(0.1),
#'       basis = list(b = list(function(u) { u }))
#'     ),
#'     list(
#'       theta = c(0.9),
#'       basis = list(b = list(function(u) { u }))
#'     )
#'   )
#' )
#'
#' # Simulate a dataset of size 100
#' set.seed(123)
#' sim_data <- simulate_data(dgp_example, 100)
#' head(sim_data)
#' 
#' @importFrom dplyr bind_rows
#' @importFrom tibble add_column
#' @importFrom stats runif rbinom
#' @export
simulate_data <- function(dgp, 
                          n, 
                          extra.details = FALSE){
  
  zvals <- sample(x = dgp$suppZ, 
                  size = n, 
                  prob = dgp$densZ, 
                  replace = TRUE)
  
  uvals <- stats::runif(n)
  
  data <- lapply(1:n,
                 function(i){
                   # D = 1 if u <= pscore(z)
                   D.temp <- uvals[i] <= dgp$find_pscore(zvals[i])
                   
                   probY.temp <- 0
                   
                   # If D == 0, use m_0 (first element of dgp$mtrs)
                   # If D == 1, use m_1 (second element of dgp$mtrs)
                   for(j in 1:length(dgp$mtrs[[D.temp + 1]]$theta)){
                     probY.temp <- probY.temp + 
                       dgp$mtrs[[D.temp + 1]]$theta[j]*dgp$mtrs[[D.temp + 1]]$basis$b[[j]](uvals[i])
                   }
                   
                   # Draw single Y, where Y ~ Bernoulli(probY)
                   Y.temp <- stats::rbinom(1, 1, probY.temp)
                   
                   return(c("y" = Y.temp,
                            "d" = D.temp,
                            "probability.y" = probY.temp))
                 }) |>
    dplyr::bind_rows() |>
    tibble::add_column(u = uvals, 
                       z = zvals)
  
  if(!extra.details){
    data$probability.y <- NULL
    data$u <- NULL
  }
  
  return(data)
  
}

#' Compute Test Statistics in Zhu (2020)
#' 
#' The function computes the test statistic for each value in \code{beta.nulls}.
#' Gurobi is used as the QP solver. Change beta.nulls to beta.null?
#' 
#' @param data A data frame or tibble that must contain the following columns:
#'   \describe{
#'     \item{\code{y}}{Binary outcome variable (0/1).}
#'     \item{\code{d}}{Binary treatment indicator (0/1).}
#'     \item{\code{z}}{Binary instrument indicator (0/1).}
#'   }
#'   
#' @param bases A list of basis functions.
#' @param beta.null A numeric vector of target parameter values for which the test 
#'  statistic is solved.  
#' @param C0 Vector of C0 values (using sample data).
#' @param C1 Matrix of C1 values (using sample data).
#' @param tau A numeric matrix representing the computed tau values.
#' @param return.gurobi Logical; if TRUE, returns raw Gurobi output.
#' @param gurobi.tries Integer; number of attempts to optimize using Gurobi.
#'  
#' @return A tibble where each row corresponds to one \code{beta.null} from the 
#'   input vector \code{beta.nulls}. The returned columns are:
#'   \itemize{
#'     \item \code{t.value}: The target parameter value used in the constraint.
#'     \item \code{test.stat}: The minimized QP objective value. 
#'       If the problem is infeasible, this is \code{NA}.
#'     \item \code{gurobi.result}: A list-column containing the raw Gurobi output
#'       for that \code{beta.null}.
#'   } 
#'   
#' @importFrom dplyr bind_rows
#' @importFrom tibble tibble
#' @importFrom gurobi gurobi
#' @importFrom cli cli_abort
#' @export   
compute_test_stat <- function(data, 
                              bases, 
                              beta.null, 
                              C0, 
                              C1, 
                              tau, 
                              return.gurobi = FALSE,
                              gurobi.tries){
  
  n <- nrow(data)  
  
  # QP formulation
  constant <- sum(C0^2)
  quad <- t(C1)%*%C1  
  linear <- matrix(-2*t(C0)%*%C1)
  
  # Lower and upper bounds for theta
  lb_theta <- rep(0, length(linear))
  ub_theta <- rep(1, length(linear))
  
  test.stat.results <- lapply(beta.null,
                              function(x){
                                
                                for(try in 1:gurobi.tries){
                                  gurobi.result <- gurobi::gurobi(list(Q=quad, 
                                                                       obj=linear, 
                                                                       modelsense="min", 
                                                                       A=t(tau), 
                                                                       rhs=x, 
                                                                       sense="=", 
                                                                       lb=lb_theta, 
                                                                       ub=ub_theta), 
                                                                  params=list(OutputFlag=0)) |>
                                    try()
                                  
                                  if(class(gurobi.result)[1] != "try-error"){
                                    break
                                  } else if(try == gurobi.tries){
                                    cli::cli_abort("Reached maximum attempts to calculate bootstrap test statistic (50) using Gurobi.")
                                  }
                                } 
                                
                                temp <- tibble::tibble("beta.null" = x)
                                
                                if(gurobi.result$status == "INFEASIBLE"){
                                  temp$test.stat <- NA
                                } else{
                                  # Add constant term to objective value
                                  gurobi.result$objval <- gurobi.result$objval + constant
                                  
                                  temp$test.stat <- n*gurobi.result$objval
                                }
                                
                                if(return.gurobi){
                                  temp$gurobi.result <- list(gurobi.result)
                                }
                                
                                return(temp)
                              }) |>
    dplyr::bind_rows()
  
  return(test.stat.results)  
}

#' Compute Bootstrap Test Statistics
#'
#' This function computes bootstrap test statistics for a given dataset using Gurobi.
#'
#' @param data A data frame or tibble containing the observed data.
#' @param beta.null Numeric vector of target parameter values.
#' @param bases List of basis functions.
#' @param C0 Vector of C0 values (using sample data).
#' @param C1 Matrix of C1 values (using sample data).
#' @param tau A numeric matrix representing the computed tau values.
#' @param gamma.lambda.df Data frame of gamma and lambda values.
#' @param number.bootstraps Number of bootstrap samples (default: 500).
#' @param bootstrap.seeds Optional vector of bootstrap seeds.
#' @param parallel Logical; whether to use parallel computation (default: TRUE).
#' @param rest.time Numeric; waiting time between optimization attempts.
#' @param dgp.in Optional; data-generating process specification.
#' @param return.gurobi Logical; if TRUE, returns raw Gurobi output.
#' @param gurobi.tries Integer; number of attempts to optimize using Gurobi.
#'
#' @return A tibble containing computed bootstrap test statistics for each combination of \code{beta.null}, \code{gamma}, and \code{lambda}.
#'
#' @importFrom parallel detectCores
#' @importFrom pbapply pblapply
#' @importFrom dplyr bind_rows
#' @export
compute_bootstrap_test_stats <- function(data,
                                         beta.null,
                                         bases,
                                         C0, 
                                         C1,
                                         tau,
                                         gamma.lambda.df,
                                         number.bootstraps = 500,
                                         bootstrap.seeds,
                                         parallel,
                                         rest.time = 0,
                                         dgp.in = NULL,
                                         return.gurobi = FALSE,
                                         gurobi.tries){
  
  
  if(is.logical(parallel)){
    if(parallel == TRUE){
      parallel <- parallel::detectCores()
    } else{
      parallel <- 1
    }
  }
  
  final.result <- pbapply::pblapply(1:number.bootstraps,
                                    function(bootstrap.number){
                                      
                                      compute_bootstrap_test_stat(bootstrap.number,
                                                                  bootstrap.seed = bootstrap.seeds[bootstrap.number],
                                                                  data,
                                                                  bases,
                                                                  beta.null,
                                                                  C0,
                                                                  C1,
                                                                  tau,
                                                                  gamma.lambda.df,
                                                                  rest.time,
                                                                  return.gurobi,
                                                                  gurobi.tries)
                                    
                                     },
                                    cl = parallel) |>
    dplyr::bind_rows()
  
  return(final.result)
}

#' Compute Bootstrap Test Statistics
#'
#' This function computes bootstrap test statistics for a given dataset using Gurobi.
#'
#' @param bootstrap.number Integer; current bootstrap number.
#' @param bootstrap.seed Integer; Used for seed to generate bootstrap data.
#' @param data A data frame or tibble containing the observed data.
#' @param bases List of basis functions.
#' @param beta.null Numeric vector of target parameter values.
#' @param C0 Vector of C0 values (using sample data).
#' @param C1 Matrix of C1 values (using sample data).
#' @param tau A numeric matrix representing the computed tau values.
#' @param gamma.lambda.df Data frame of gamma and lambda values.
#' @param rest.time Numeric; waiting time between optimization attempts.
#' @param return.gurobi Logical; if TRUE, returns raw Gurobi output.
#' @param gurobi.tries Integer; number of attempts to optimize using Gurobi.
#'
#' @return A tibble containing computed bootstrap test statistics for each combination of \code{beta.null}, \code{gamma}, and \code{lambda}.
#'
#' @importFrom dplyr bind_rows
#' @importFrom gurobi gurobi
#' @importFrom tibble tibble
#' @importFrom cli cli_abort
#' @export
compute_bootstrap_test_stat <- function(bootstrap.number,
                                        bootstrap.seed,
                                        data,
                                        bases,
                                        beta.null,
                                        C0,
                                        C1,
                                        tau,
                                        gamma.lambda.df,
                                        rest.time,
                                        return.gurobi,
                                        gurobi.tries){
  
  n <- nrow(data)
  
  set.seed(bootstrap.seed)
  
  # If bootstrap == TRUE, then we draw with replacement and 
  # uniform probability from original sample B
  bootstrap.indices <- sample(1:n, n, replace = TRUE)
  bootstrap.data <- data[bootstrap.indices,]
  
  # Compute C0.star and C1.star (Reference's Page 6)
  C0.C1.star.list <- compute_C0_C1(bootstrap.data, bases)
  C0.star <- C0.C1.star.list$C0
  C1.star <- C0.C1.star.list$C1
  
  result <- lapply(beta.null,
                   function(x){
                     apply(gamma.lambda.df,
                           MARGIN = 1,
                           function(y){
                             gamma_n <- y["gamma"]
                             lambda_n <- y["lambda"]
                             
                             # Compute quadratic term matrix and linear term vector
                             # Note: we don't need the constant term here because we are interested in arg min 
                             # over \theta. Adding constant to objective function does not affect arg min.
                             QLC <- compute_QLC_nQstar_n(n, 
                                                         C0, 
                                                         C0.star, 
                                                         C1, 
                                                         C1.star, 
                                                         gamma_n, 
                                                         lambda_n)
                             quad <- QLC$quadratic 
                             linear <- QLC$linear
                             
                             # Solve optimization problem
                             for(try in 1:gurobi.tries){
                               gurobi.result <- gurobi::gurobi(list(Q=quad, 
                                                                    obj=linear, 
                                                                    modelsense="min", 
                                                                    A=t(tau), 
                                                                    rhs=x, 
                                                                    sense="=", 
                                                                    lb=rep(0, length(C0)), # Declare lower and upper bounds for theta for bernstein basis 
                                                                    ub=rep(1, length(C0))), 
                                                               params=list(OutputFlag=0,
                                                                           threads = 1)) |>
                                 try()
                               
                               if(class(gurobi.result)[1] != "try-error"){
                                 break
                               } else if(try == 50){
                                 cli::cli_abort("Reached maximum attempts to calculate bootstrap test statistic (50) using Gurobi.")
                               }
                               
                               Sys.sleep(rest.time)
                             }
                             
                             temp <- tibble::tibble(bootstrap.number,
                                                    bootstrap.seed,
                                                    beta.null = x,
                                                    gamma = gamma_n,
                                                    lambda = lambda_n)
                             
                             if(gurobi.result$status == "INFEASIBLE"){
                               temp$bootstrap.test.stat <- NA
                             } else{
                               solution.theta <- gurobi.result$x
                               
                               # Compute new quadratic term matrix, linear term vector, and constant term
                               # This is for Reference's Eqn. 10
                               QLC.new <- compute_QLC_nQstar_n(n, C0, C0.star, C1, C1.star, gamma_n, 0)
                               quad.new <- QLC.new$quadratic 
                               linear.new <- QLC.new$linear
                               constant.new <- QLC.new$constant
                               
                               TStar.n <- as.numeric(t(solution.theta)%*%quad.new%*%solution.theta + 
                                                       linear.new%*%solution.theta + 
                                                       constant.new)
                               
                               # Add contender bootstrap test statistic to TStar_n 
                               # No need to multiply TStar_n by n here since already did in compute_QLC_nQstar_n().
                               temp$bootstrap.test.stat <- TStar.n
                             }
                             
                             if(return.gurobi){
                               temp$gurobi.result <- list(gurobi.result)
                             }
                             
                             return(temp)  
                           },
                           simplify = FALSE) |>
                       dplyr::bind_rows()
                   }) |>
    dplyr::bind_rows()
  
  return(result)
  
}

#' Compute C0 and C1 Matrices
#'
#' Computes the C0 and C1 matrices as defined in Torgovitsky (2020).
#'
#' @param data A data frame containing outcome \code{y}, treatment \code{d}, and covariate \code{z}.
#' @param bases A list of basis functions.
#'
#' @return A list with elements:
#'   \itemize{
#'     \item \code{C0}: Vector of estimated C0 values.
#'     \item \code{C1}: Matrix of estimated C1 values.
#'   }
#'
#' @importFrom dplyr group_by summarise %>%
#' @importFrom cli cli_abort
#' @importFrom rlang .data
compute_C0_C1 <- function(data, bases){
  
  if(!all(c("y", "d", "z") %in% colnames(data))){
    cli::cli_abort("`data` does not include columns 'y', 'd', and 'z'.")
  }
  
  pscoreZ <- data |>
    dplyr::group_by(.data$z) |>
    dplyr::summarise(propScore = mean(.data$d), .groups = "drop")
  
  # Build B function for each observations in data
  B.matrix <- apply(data.frame(data),
                    MARGIN = 1,
                    function(row){
                   
                      d.i <- row["d"]
                      z.i <- row["z"]
                      
                      pscoreZi <- pscoreZ$propScore[pscoreZ$z == z.i]
                   
                      return(compute_B(z.i, 
                                       d.i, 
                                       pscoreZi,
                                       bases))
                    })
  
  C0 <- B.matrix %*% matrix(data$y)
  
  C1 <- Reduce(`+`,
               apply(B.matrix,
                     MARGIN = 2,
                     function(x){
                       outer(x, x)
                     },
                     simplify = FALSE))
  
  return(list(C0 = C0/nrow(B.matrix), 
              C1 = C1/nrow(B.matrix)))
}

#' Compute Quadratic, Linear, and Constant Terms
#'
#' Computes the QLC terms used in quadratic optimization for hypothesis testing.
#'
#' @param n Sample size.
#' @param C0 Vector of C0 values (using sample data).
#' @param C0.star Vector of C0 values (using bootstrap data).
#' @param C1 Matrix of C1 values (using sample data).
#' @param C1.star Matrix of C1 values (using bootstrap data).
#' @param gamma_n Scalar gamma value.
#' @param lambda_n Scalar lambda value.
#'
#' @return A list containing:
#'   \itemize{
#'     \item \code{quadratic}: Quadratic term matrix.
#'     \item \code{linear}: Linear term vector.
#'     \item \code{constant}: Constant term scalar.
#'   }
#'   
compute_QLC_nQstar_n = function(n, C0, C0.star, C1, C1.star, gamma_n, lambda_n){
  
  A <- C0.star + (gamma_n - 1)*C0
  B <- C1.star + (gamma_n - 1)*C1
  
  # Compute quadratic term matrix
  quadratic <- n*t(B)%*%B + lambda_n*t(C1)%*%C1
  
  # Compute linear term vector
  linear <- -2*(n*t(A)%*%(B) + lambda_n*t(C0)%*%C1)
  
  # Compute constant term
  constant <- n*sum(A^2) + lambda_n*sum(C0^2)
  
  return(list(quadratic = quadratic, linear = linear, constant = constant))
}

#' Compute Tau Matrix
#'
#' This function calculates the tau matrix as discussed in Torgovitsky (2023).
#' It determines weighted averages of basis functions over the sample or a given DGP.
#'
#' @param data A data frame containing the observed sample data. Should include columns:
#'   \itemize{
#'     \item \code{y}: Outcome variable.
#'     \item \code{d}: Treatment indicator (0 or 1).
#'     \item \code{z}: Covariate.
#'   }
#'   Set to NULL if using \code{dgp}.
#' @param dgp A data-generating process object containing the probability distribution of \code{Z}.
#'   Set to NULL if using \code{data}.
#' @param bases A list of basis functions, where:
#'   \itemize{
#'     \item \code{bases[[1]]$ib}: Basis functions for the control group (\code{d = 0}).
#'     \item \code{bases[[2]]$ib}: Basis functions for the treated group (\code{d = 1}).
#'   }
#' @param target.parameter The treatment effect parameter to estimate (e.g., ATT, ATE, LATE).
#' @param late.lb Lower bound for LATE estimation (required when target.parameter = "LATE")
#' @param late.ub Upper bound for LATE estimation (required when target.parameter = "LATE")
#'
#' @return A numeric matrix representing the computed tau values.
#'
#' @importFrom ivprte compute_average_weights
compute_tau <- function(data, 
                        dgp, 
                        bases, 
                        target.parameter, 
                        late.lb, 
                        late.ub){
  
  if((is.null(data) & is.null(dgp)) | (!is.null(data) & !is.null(dgp))){
    stop("The `compute_tau()` function accepts sample data or dgp (exclusive) as parameters.")
  }
  
  if(!is.null(data)){
    weights.df <- ivprte::compute_average_weights(target.parameter, data = data, dgp = NULL, late.lb, late.ub)
  } else{
    weights.df <- ivprte::compute_average_weights(target.parameter, data = NULL, dgp = dgp, late.lb, late.ub)
  }
  
  tau.d0 <- sapply(bases[[1]]$ib,
                   function(integrate){
                     apply(weights.df,
                           MARGIN = 1,
                           function(x){
                             
                             avg.weight.d0 <- x["avgWeightD0"]
                             
                             if(avg.weight.d0 == 0){
                               return(0)
                             }
                             
                             u.start <- x["uStart"]
                             u.end <- x["uEnd"]
                             
                             return(avg.weight.d0*integrate(u.start, u.end))
                             
                           }) |>
                       sum()
                   })
  
  tau.d1 <- sapply(bases[[2]]$ib,
                   function(integrate){
                     apply(weights.df,
                           MARGIN = 1,
                           function(x){
                             
                             avg.weight.d1 <- x["avgWeightD1"]
                             
                             if(avg.weight.d1 == 0){
                               return(0)
                             }
                             
                             u.start <- x["uStart"]
                             u.end <- x["uEnd"]
                             
                             return(avg.weight.d1*integrate(u.start, u.end))
                             
                           }) |>
                       sum()
                   })
  
  return(matrix(c(tau.d0, tau.d1), ncol = 1))
}

#' Compute B Matrix for Basis Functions
#'
#' This function computes the B matrix based on the treatment indicator \code{d},
#' the propensity score \code{pscoreZ}, and a set of basis functions.
#'
#' @param z The covariate value (unused in the function but included for consistency).
#' @param d Binary treatment indicator (0 or 1).
#' @param pscoreZ The propensity score for treatment assignment given \code{Z}.
#' @param bases A list of basis functions, where:
#'   \itemize{
#'     \item \code{bases[[1]]$ib}: Basis functions for the control group (\code{d = 0}).
#'     \item \code{bases[[2]]$ib}: Basis functions for the treated group (\code{d = 1}).
#'   }
#'
#' @return A numeric vector representing the computed B matrix for the given inputs.
#' 
compute_B = function(z, d, pscoreZ, bases){
  
  if(d == 0){
    
    B <- c(sapply(bases[[1]]$ib,
                  function(x){
                    x(pscoreZ, 1)/(1-pscoreZ)
                  }),
           rep(0, length(bases[[2]]$ib))
           )
    
  } else{
    
    B <- c(rep(0, length(bases[[1]]$ib)),
           sapply(bases[[2]]$ib,
                    function(x){
                      x(0, pscoreZ)/pscoreZ
                      })
           )
    
  }
  
  return(B)
}

#' Generate Gamma-Lambda Data Frame
#'
#' Generates a data frame of gamma and lambda values for use in optimization.
#'
#' @param kappas Vector of kappa values.
#' @param sample.size Integer; sample size.
#' @param MB Integer; 1 for endpoint method, 2 for equally spaced gamma values.
#'
#' @return A data frame with columns \code{gamma}, \code{lambda}, and \code{kappa}.
#'
#' @importFrom dplyr mutate select
#' @importFrom rlang .data
#' @export
compute_gamma_lambda_df <- function(kappas, sample.size, MB = 1){
  
  # For each kappa value, we generate a data.frame() with the corresponding gamma and lambda values
  results <- lapply(kappas,
                    function(kappa){
                      # At this point, we know that MB must be 1 or 2.
                      if(MB == 1){ 
                        temp <- data.frame(gamma = c(sqrt(kappa/sample.size), 0),
                                           lambda = c(0, kappa),
                                           kappa)
                      } else{
                        temp <- data.frame(gamma = seq(sqrt(kappa/sample.size), 0, -sqrt(kappa/sample.size)/30),
                                           kappa) |>
                          dplyr::mutate(lambda = round(kappa - sample.size*(.data$gamma^2), 10)) |>
                          dplyr::select(.data$gamma, 
                                        .data$lambda, 
                                        .data$kappa)
                      }
                      
                      return(temp)
                    })
  
  return(unique(do.call('rbind', results)))
}

#' Shape-Restricted Inference for Treatment Effects
#'
#' This function performs shape-restricted inference for a given target parameter.
#'
#' @param data A data frame containing columns `y` (outcome), `d` (treatment), and `z` (instrument).
#' @param target.parameter The treatment effect parameter to estimate (e.g., ATT, ATE, LATE).
#' @param alpha Significance level for confidence intervals (default: 0.05).
#' @param bases A list of basis functions used in the estimation.
#' @param beta.null A vector of null values for the hypothesis test.
#' @param number.bootstraps Number of bootstrap samples (default: 500).
#' @param bootstrap.seeds Optional vector of bootstrap seeds.
#' @param gamma.lambda.df Optional data frame of gamma-lambda values.
#' @param kappa Regularization parameter (if `gamma.lambda.df` is not provided).
#' @param MB Bootstrap method indicator (default: 1).
#' @param dgp.in Optional; data-generating process specification. If included, computes population \eqn{\tau} rather than estimated \eqn{\tau}.
#' @param parallel Logical; whether to use parallel computation (default: TRUE).
#' @param return.gurobi Logical; whether to return raw Gurobi output (default: FALSE).
#' @param gurobi.tries Number of attempts for Gurobi solver (default: 50).
#' @param late.lb Lower bound for LATE estimation (required when target.parameter = "LATE")
#' @param late.ub Upper bound for LATE estimation (required when target.parameter = "LATE")
#'
#' @return An object of class `shapeinf` containing:
#'   \itemize{
#'     \item \code{ci.bounds}: Confidence interval for the target parameter.
#'     \item \code{beta.null.test.detailed}: Detailed test statistics for `beta.null` values.
#'     \item \code{bootstrap.detailed}: Detailed bootstrap test statistics.
#'     \item \code{data}: The input data.
#'     \item \code{bases}: The basis functions used.
#'     \item \code{tau}: Computed tau values.
#'     \item \code{C0}, \code{C1}: Computed C0 and C1 matrices.
#'   }
#'
#' @importFrom cli cli_abort cli_alert_info
#' @importFrom dplyr left_join mutate group_by summarise
#' @importFrom rlang .data
#' @importFrom Matrix Matrix
#' @export
shapeinf <- function(data,
                     target.parameter,
                     alpha = 0.05,
                     bases,
                     beta.null,
                     number.bootstraps = 500,
                     bootstrap.seeds = NULL,
                     gamma.lambda.df = NULL,
                     kappa = NULL,
                     MB = 1,
                     dgp.in = NULL,
                     parallel = TRUE,
                     return.gurobi = FALSE,
                     gurobi.tries = 50,
                     late.lb = NULL,
                     late.ub = NULL){
  
  
  ### Random useless Matrix package code.
  ### Used to resolve devtools::check() error/note
  
  junk <- Matrix::Matrix(1)
  
  call <- match.call()
  call$target.parameter <- eval(call$target.parameter, envir = parent.frame())
  call$alpha <- eval(call$alpha, envir = parent.frame())
  
  if(!is.null(late.lb)){
    call$late.lb <- eval(call$late.lb, envir = parent.frame())
  }
  
  if(!is.null(late.ub)){
    call$late.ub <- eval(call$late.ub, envir = parent.frame())
  }
  
  # Future warm start option
  # if(!is.null(warm.start)){
  #   if(class(warm.start) != "shapeinf"){
  #     cli::cli_abort("The warm.start parameter only accepts objects of class 'shapeinf'.")
  #   }
  #   
  #   og.target.parameter <- rownames(warm.start$ci.bounds)
  #   og.bases <- warm.start$bases
  #   
  #   if(og.target.parameter != target.parameter){
  #     cli::cli_abort("The target parameter of the warm start object must be equal to the current target parameter.")
  #   }
  #   
  #   if(!all.equal(og.bases, bases)){
  #     cli::cli_abort("The bases parameter of the warm start object must be equal to the current bases parameter.")
  #   }
  #   
  #   if(og.target.parameter == "LATE"){
  #     og.late.lb <- warm.start$late.lb
  #     og.late.ub <- warm.start$late.ub
  #     
  #     if(og.late.lb != late.lb || og.late.ub != late.ub){
  #       cli::cli_abort("The late.lb and late.ub parameters of the warm start object must be equal to the current late.lb and late.ub parameters.")
  #     }
  #   }
  #   
  #   if(!is.null(data) && !all.equal(data[, c('y', 'd', 'z')], warm.start$data[, c('y', 'd', 'z')])){
  #     cli::cli_abort("The data parameter of the warm start object does not equal the current data.")
  #   }
  # }
  
  sample.size <- nrow(data)
  
  if(is.null(gamma.lambda.df)){
    
    if(MB != 1 & MB != 2){
      cli::cli_abort("`MB` must be set to '1' or '2'.")
    }
    
    if(is.null(kappa)){
      kappa <- sample.size/log(sample.size)
    }
    
    gamma.lambda.df <- compute_gamma_lambda_df(kappas = kappa, 
                                               sample.size,
                                               MB = MB)
    
  } else if(!((is.data.frame(gamma.lambda.df) || inherits(gamma.lambda.df, "tbl_df")) &&
            all(c("gamma", "lambda") %in% colnames(gamma.lambda.df)))){ 
    
    # Check if gamma.lambda.df is data.frame or tibble with columns gamma and lambda
    cli::cli_abort("gamma.lambda.df must be a data.frame or tibble with column names that include 'gamma' and 'lambda'.")
    
  } 
  
  # Check if data is data.frame and contains columns y, d, and z
  if(!((is.data.frame(data) || inherits(data, "tbl_df")) && 
       all(c("y", "d", "z") %in% colnames(data)))){
    cli::cli_abort("data must be a data.frame or tibble with column names that include 'y', 'd', and 'z'.")
  }
  
  if(!is.null(bootstrap.seeds) && 
     length(bootstrap.seeds) != number.bootstraps){
    cli::cli_abort("If seeds for each bootstrap are used, the number of bootstrap seeds must be equal to the number of bootstraps.")
  } else if(is.null(bootstrap.seeds)){
    bootstrap.seeds <- sample.int(1e8 + number.bootstraps, 
                                  size = number.bootstraps)[1:number.bootstraps]
  }
  
  # Compute C0 and C1
  C0.C1.list <- compute_C0_C1(data, bases)
  
  C0 <- C0.C1.list$C0
  C1 <- C0.C1.list$C1
  
  if(is.null(dgp.in)){
    tau <- compute_tau(target.parameter, 
                       data = data, 
                       dgp = NULL, 
                       bases = bases, 
                       late.lb, 
                       late.ub)
  } else{
    tau <- compute_tau(target.parameter, 
                       data = NULL,
                       dgp = dgp.in, 
                       bases = bases, 
                       late.lb, 
                       late.ub)
  }
  
  cli::cli_alert_info("Computing test statistics...")
  test.stats <- compute_test_stat(data, 
                                  bases, 
                                  beta.null,
                                  C0,
                                  C1,
                                  tau,
                                  return.gurobi = return.gurobi,
                                  gurobi.tries = gurobi.tries)
  
  cli::cli_alert_info("Computing bootstrap test statistics...")
  bootstrap.test.stats <- compute_bootstrap_test_stats(data,
                                                       beta.null,
                                                       bases,
                                                       C0,
                                                       C1,
                                                       tau,
                                                       gamma.lambda.df,
                                                       number.bootstraps,
                                                       bootstrap.seeds = bootstrap.seeds,
                                                       parallel = parallel,
                                                       return.gurobi = return.gurobi,
                                                       gurobi.tries = gurobi.tries)
  
  cli::cli_alert_info("Computing the confidence interval through test inversion...")
  
  p.values <- test.stats |>
    dplyr::left_join(bootstrap.test.stats,
                     by = "beta.null") |>
    dplyr::mutate(bootstrap.test.stat.larger = .data$bootstrap.test.stat >= .data$test.stat) |>
    dplyr::group_by(.data$beta.null) |>
    dplyr::summarise(p.value = mean(.data$bootstrap.test.stat.larger))
  
  beta.null.test.detailed <- test.stats |>
    dplyr::left_join(p.values,
                     by = "beta.null")
  
  lower.bound.index <- which(beta.null.test.detailed$p.value >= alpha) |>
    min() |>
    suppressWarnings()
  
  # Protect against scenario where we need 0 index or no p.value is >= than alpha
  if(lower.bound.index == 1 || lower.bound.index == Inf){
    cli::cli_alert_warning("Unable to locate a lower bound. Try including smaller values in beta.null.")
    lower.bound <- NA
  } else{
    slope <- (beta.null.test.detailed$p.value[lower.bound.index] - beta.null.test.detailed$p.value[lower.bound.index - 1])/
      (beta.null.test.detailed$beta.null[lower.bound.index] - beta.null.test.detailed$beta.null[lower.bound.index - 1])
    
    lower.bound <- (alpha - beta.null.test.detailed$p.value[lower.bound.index])/slope + 
      beta.null.test.detailed$beta.null[lower.bound.index]
  }
  
  upper.bound.index <- which(beta.null.test.detailed$p.value >= alpha) |>
    max() |>
    suppressWarnings()
  
  # Protect against scenario where we need sample.size + 1 index or no p.value is >= than alpha
  if(upper.bound.index == sample.size || upper.bound.index == Inf){
    cli::cli_alert_warning("Unable to locate an upper bound. Try including larger values in beta.null.")
    upper.bound <- NA
  } else{
    slope <- (beta.null.test.detailed$p.value[upper.bound.index] - beta.null.test.detailed$p.value[upper.bound.index + 1])/
      (beta.null.test.detailed$beta.null[upper.bound.index] - beta.null.test.detailed$beta.null[upper.bound.index + 1])
    
    upper.bound <- (alpha - beta.null.test.detailed$p.value[upper.bound.index])/slope + 
      beta.null.test.detailed$beta.null[upper.bound.index]
  }
  
  ci.bounds <- matrix(data = c(lower.bound, upper.bound),
                      nrow = 1,
                      ncol = 2,
                      dimnames = list(target.parameter,
                                      c("lower bound", "upper bound")))
  
  result <- structure(list(call = call,
                           ci.bounds = ci.bounds,
                           beta.null.test.detailed = beta.null.test.detailed,
                           bootstrap.detailed = bootstrap.test.stats,
                           data = data,
                           bases = bases,
                           tau = tau,
                           C0 = C0,
                           C1 = C1),
                      class = "shapeinf")
  
  return(result)
}

#' Expected B B' Matrix Calculation
#'
#' Computes the expected outer product of B matrices given a DGP and basis functions.
#'
#' @param dgp A list describing the data-generating process.
#' @param bases A list of basis functions.
#' @return A matrix representing the expected outer product of B matrices.
#' @export
expectBBPrime = function(dgp, bases){
  
  dim <- length(bases[[1]]$ib) + length(bases[[2]]$ib)
  
  # Initialize result matrix
  result <- matrix(0, ncol = dim, nrow = dim)
  
  for(z in dgp$suppZ){
    
    probZ <- dgp$find_density(z)
    
    for(d in 0:1){
      pscoreZ <- dgp$find_pscore(z)
      probDConditionalOnZ <- ifelse(d == 1, pscoreZ, 1-pscoreZ)
      
      result <- result + 
        probZ*probDConditionalOnZ*outer(compute_B(z, d, pscoreZ, bases), 
                                        compute_B(z, d, pscoreZ, bases)) 
    }
  }
  return(result)
}

#' Expected Outcome Given D and Z
#'
#' Computes the expected outcome \deqn{E[Y \mid D, Z]}{E[Y | D, Z]} using MTR coefficients.
#'
#' @param z The instrument value.
#' @param d The treatment indicator (0 or 1).
#' @param dgp A list describing the data-generating process.
#' @return Expected outcome given D and Z.
#' @export
expectYConditionalOnDZ = function(z, d, dgp){
  
  theta <- c(dgp$mtrs[[1]]$theta, dgp$mtrs[[2]]$theta)
  bases <- list(dgp$mtrs[[1]]$basis, dgp$mtrs[[2]]$basis)
  
  return(sum(theta*compute_B(z, d, dgp$find_pscore(z), bases)))
}

#' Expected B Y Product Calculation
#'
#' Computes the expected product of B and Y given a DGP and basis functions.
#'
#' @param dgp A list describing the data-generating process.
#' @param bases A list of basis functions.
#' @return A vector representing the expected B Y product.
#' @export
expectBY = function(dgp, bases){
  
  length.bases <- length(bases[[1]]$ib) + length(bases[[2]]$ib)
  
  result <- rep(0, length.bases)
  
  for(z in dgp$suppZ){
    for(d in 0:1){
      probD <- ifelse(d == 1, dgp$find_pscore(z), 1 - dgp$find_pscore(z))
      probZ <- dgp$find_density(z)
      
      result <- result + probZ * probD * compute_B(z, d, dgp$find_pscore(z), bases) * 
        expectYConditionalOnDZ(z, d, dgp)
    }
  }
  return(result)
}

#' Print Method for `shapeinf` Objects
#'
#' @param x An object of class `shapeinf`.
#' @param ci.digits Number of digits to display in confidence intervals.
#' @param ... Additional arguments (not currently used).
#' 
#' @importFrom cli style_bold col_blue
#' @export
print.shapeinf <- function(x, ci.digits = 4, ...){
  
  cat(cli::style_bold(cli::col_blue("\nZhu Shape-Restricted Test\n")))
  
  cat(cli::style_bold(cli::col_blue("\nCall:")),
      paste(deparse(x$call), sep = "\n", collapse = "\n"), 
      "\n", 
      sep = "")
  
  cat(cli::style_bold(cli::col_blue("\nTarget Parameter: ")),
      switch(x$call$target.parameter,
             "AUO" = "Average Untreated Outcome (AUO)",
             "ATO" = "Average Treated Outcome (ATO)",
             "ATE" = "Average Treatment Effect (ATE)",
             "ATT" = "Average Treatment on the Treated (ATT)",
             "ATU" = "Average Treatment on the Untreated (ATU)",
             "LATE" = sprintf("Local Average Treatment Effect for U in (%s, %s] (LATE(%s, %s))", 
                              format(x$call$late.lb), 
                              format(x$call$late.ub), 
                              format(x$call$late.lb), 
                              format(x$call$late.ub)),
             x$call$target.parameter),
      "\n\n",
      sep = "")
  
  cat(cli::style_bold(cli::col_blue(100 * (1-x$call$alpha), "% Confidence Interval:")),
      paste0("(", signif(x$ci.bounds[1], ci.digits), ", ", signif(x$ci.bounds[2], ci.digits), ")"))
}


#' Summary Method for `shapeinf` Objects
#'
#' @param object An object of class `shapeinf`.
#' @param alpha Optional new significance level.
#' @param ... Additional arguments (not currently used).
#' 
#' @importFrom cli cli_alert_warning
#' @export
summary.shapeinf <- function(object, alpha = NULL, ...){
  
  if(!is.null(alpha) && alpha != object$call$alpha){
    object$call$alpha <- alpha
    
    sample.size <- nrow(object$data)
    
    lower.bound.index <- which(object$beta.null.test.detailed$p.value >= alpha) |>
      min() |>
      suppressWarnings()
    
    # Protect against scenario where we need 0 index or no p.value is >= than alpha
    if(lower.bound.index == 1 || lower.bound.index == Inf){
      cli::cli_alert_warning("Unable to locate a lower bound. Try including smaller values in beta.null.")
      lower.bound <- NA
    } else{
      slope <- (object$beta.null.test.detailed$p.value[lower.bound.index] - object$beta.null.test.detailed$p.value[lower.bound.index - 1])/
        (object$beta.null.test.detailed$beta.null[lower.bound.index] - object$beta.null.test.detailed$beta.null[lower.bound.index - 1])
      
      lower.bound <- (alpha - object$beta.null.test.detailed$p.value[lower.bound.index])/slope + 
        object$beta.null.test.detailed$beta.null[lower.bound.index]
    }
    
    upper.bound.index <- which(object$beta.null.test.detailed$p.value >= alpha) |>
      max() |>
      suppressWarnings()
    
    # Protect against scenario where we need sample.size + 1 index or no p.value is >= than alpha
    if(upper.bound.index == sample.size || upper.bound.index == Inf){
      cli::cli_alert_warning("Unable to locate an upper bound. Try including larger values in beta.null.")
      upper.bound <- NA
    } else{
      slope <- (object$beta.null.test.detailed$p.value[upper.bound.index] - object$beta.null.test.detailed$p.value[upper.bound.index + 1])/
        (object$beta.null.test.detailed$beta.null[upper.bound.index] - object$beta.null.test.detailed$beta.null[upper.bound.index + 1])
      
      upper.bound <- (alpha - object$beta.null.test.detailed$p.value[upper.bound.index])/slope + 
        object$beta.null.test.detailed$beta.null[upper.bound.index]
    }
    
    ci.bounds <- matrix(data = c(lower.bound, upper.bound),
                        nrow = 1,
                        ncol = 2,
                        dimnames = list(object$call$target.parameter,
                                        c("lower bound", "upper bound")))
    
    object$ci.bounds <- ci.bounds
  }
  
  print(object)  
}




