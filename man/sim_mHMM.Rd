% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sim_mHMM.R
\name{sim_mHMM}
\alias{sim_mHMM}
\title{Simulate data using a multilevel hidden Markov model}
\usage{
sim_mHMM(n_t, n, m, q_emiss, gamma, emiss_distr, beta = NULL,
  xx_vec = NULL, var_gamma = 1, var_emiss = 1,
  return_ind_par = FALSE)
}
\arguments{
\item{n_t}{The length of the observed sequence to be simulated for each
subject. To only simulate subject specific transition probability matrices
gamma and emission distributions (and no data), set \code{t} to 0.}

\item{n}{The number of subjects for which data is simulated.}

\item{m}{The number of hidden states in the HMM used for
simulating data.}

\item{q_emiss}{The number of categories of the simulated observations.}

\item{gamma}{A matrix with \code{m} rows and \code{m} columns containing the
average population transition probability matrix used for simulating the
data. That is, the probability to switch from hidden state \emph{i} (row
\emph{i}) to hidden state \emph{j} (column  \emph{j}).}

\item{emiss_distr}{A matrix with \code{m} rows and \code{q_emiss} columns
containing the average population emission distribution of the
(categorical) observations given the hidden states. That is, the
probability of observing category \emph{k} (column \emph{k}) in state
\emph{i} (row \emph{i}).}

\item{beta}{List of two matrices containing the regression parameters to
predict \code{gamma} and/or \code{emiss_distr} in combination with
\code{xx_vec} using multinomial logistic regression. The first matrix is
used to predict the transition probability matrix \code{gamma}. The
second matrix is used to predict the emission distribution
\code{emiss_distr} of the dependent variable. In both matrices, one
regression parameter is specified for each element in \code{gamma} and
\code{emiss_distr}, with the following exception. The first element in each
row of \code{gamma} and/or \code{emiss_distr} is used as reference category
in the multinomial logistic regression. As such, no regression parameters
can be specified for these parameters. Hence, the first matrix in the list
\code{beta} to predict \code{gamma} consist of a matrix with the number of
rows equal to \code{m} and the number of columns equal to \code{m} - 1. The
second matrix in the list \code{beta} to predict \code{emiss_distr} consist
of a matrix with the number of rows equal to \code{m} and the number of
columns equal to \code{q_emiss} - 1. See \emph{details} for more information.
Note that if \code{beta} is specified, \code{xx_vec} has to be specified as
well. If \code{beta} is omitted completely, \code{beta} defaults to NULL,
resembling no prediction of \code{gamma} or \code{emiss_distr} using
covariates. One of the two elements in the list can also be left empty
(i.e., set to \code{NULL}) to signify that either the transition
probability matrix or a specific emission distribution is not predicted by
covariates.}

\item{xx_vec}{List of two vectors containing the covariate(s) to predict
\code{gamma} and/or \code{emiss_distr} using the regression parameters
specified in \code{beta}. The covariate used to predict \code{gamma} and
\code{emiss_distr} can either be the same covariate, two different
covariates, or a covariate for one element and none for the other. At this
point, it is only possible to use one covariate for both \code{gamma} and
\code{emiss_distr}. The first vector of the list \code{xx_vec} is used to
predict the transition matrix. The second vector of the list \code{xx_vec}
is used to predict the emission distribution of the dependent variable. For
both vectors, the number of observations should be  equal to \code{n} the
number of subjects to be simulated. If \code{xx_vec} is omitted completely,
\code{xx_vec} defaults to NULL, resembling no covariates at all. One of the
two elements in the list can also be left empty (i.e., set to \code{NULL})
to signify that either the transition probability matrix or the emission
distribution is not predicted by covariates.}

\item{var_gamma}{An integer denoting the variance between subjects in the
transition probability matrix. Note that this value corresponds to the
variance of the parameters of the multinomial distribution (i.e., the
intercepts of the regression equation of the multinomial distribution used
to sample the transition probability matrix), see details below. In
addition, only one variance value can be specified for the complete
transition probability matrix, hence the variance is assumed fixed across
all components. The default equals 1, which corresponds to quite some
variation between subjects. A less extreme value would be 0.5. If one wants
to simulate data from exactly the same HMM for all subjects, var_gamma
should be set to 0.}

\item{var_emiss}{An integer denoting the variance between subjects in the
emission distribution. Note that this value corresponds to the variance of
the parameters of the multinomial distribution (i.e., the intercepts of the
regression equation of the multinomial distribution used to sample the
components of the emission distribution), see details below.  In addition,
only one variance value can be specified for the complete emission
distribution, hence the variance is assumed fixed across all components.
The default equals 1, which corresponds to quite some variation between
subjects. A less extreme value would be 0.5. If one wants to simulate data
from exactly the same HMM for all subjects, var_emiss should be set to
0.}

\item{return_ind_par}{A logical scalar. Should the subject specific
transition probability matrix \code{gamma} and emission probability matrix
\code{emiss_distr} be returned by the function (\code{return_ind_par =
TRUE}) or not (\code{return_ind_par = FALSE}). The default equals
\code{return_ind_par = FALSE}.}
}
\value{
The following components are returned by the function \code{sim_mHMM}:
\describe{
  \item{\code{states}}{A matrix containing the simulated hidden state
  sequences, with one row per hidden state per subject. The first column
  indicates subject id number. The second column contains the simulated
  hidden state sequence, consecutively for all subjects. Hence, the id number
  is repeated over the rows (with the number of repeats equal to the length
  of the simulated hidden state sequence \code{T} for each subject).}
  \item{\code{obs}}{A matrix containing the simulated observed outputs, with
  one row per simulated observation per subject. The first column indicates
  subject id number. The second column contains the simulated observation
  sequence, consecutively for all subjects. Hence, the id number is repeated
  over rows (with the number of repeats equal to the length of the simulated
  observation sequence \code{T} for each subject).}
  \item{\code{gamma}}{A list containing \code{n} elements with the simulated
  subject specific transition probability matrices \code{gamma}. Only
  returned if \code{return_ind_par} is set to \code{TRUE}.}
  \item{\code{emiss_distr}}{A list containing \code{n} elements with the
  simulated subject specific emission probability matrices
  \code{emiss_distr}. Only returned if \code{return_ind_par} is set to
  \code{TRUE}.}
}
}
\description{
\code{sim_mHMM} simulates data for multiple subjects, for which the data have
categorical observations that follow a hidden Markov model (HMM) with an
multilevel structure. The multilevel structure implies that each subject is
allowed to have its own set of parameters, and that the parameters at the
subject level (level 1) are tied together by a population distribution at
level 2 for each of the corresponding parameters. The shape of the population
distribution for each of the parameters is a normal (i.e., Gaussian)
distribution. In addition to (natural and/or unexplained) heterogeneity
between subjects, the subjects parameters can also depend on a (set of)
covariate(s).
}
\details{
In simulating the data, having a multilevel structure means that the
parameters for each subject are sampled from the population level
distribution of the corresponding parameter. The user specifies the
population distribution for each parameter: the average population transition
probability matrix and its variance, and the average population emission
distribution and its variance. For now, the variance is assumed fixed for all
components of the transition probability matrix and for all components of the
emission distribution, and the simulated data can only consist of one
dependent variable. In addition, at this point only one dependent variable can
be simulated. That is, the hidden Markov model is a univariate hidden Markov
model.

Note: the subject specific) initial state distributions (i.e., the
probability of each of the states at the first time point) needed to simulate
the data are obtained from the stationary distributions of the subject
specific transition probability matrices gamma.

\code{beta}: As the first element in each row of \code{gamma} is used as reference category
  in the multinomial logistic regression, the first matrix in the list \code{beta} used to predict
transition probability matrix \code{gamma} has a number of rows equal to \code{m}
and the number of columns equal to \code{m} - 1. The first element in the
first row corresponds to the probability of switching from state one to state
two. The second element in the first row corresponds to the probability of
switching from state one to state three, and so on. The last element in the
first row corresponds to the probability of switching from state one to the
last state. The same principle holds for the second matrix in the list
\code{beta} used to predict the emission distribution \code{emiss_distr}: the first element
in the first row corresponds to the probability of observing category two in
state one. The second element in the first row corresponds to the probability
of observing category three is state one, and so on. The last element in the
first row corresponds to the probability of observing the last category in
state one.
}
\examples{
# simulating data for 10 subjects with each 100 observations
n_t     <- 100
n       <- 10
m       <- 3
q_emiss <- 4
gamma   <- matrix(c(0.8, 0.1, 0.1,
                  0.2, 0.7, 0.1,
                  0.2, 0.2, 0.6), ncol = m, byrow = TRUE)
emiss_distr <- matrix(c(0.5, 0.5, 0.0, 0.0,
                        0.1, 0.1, 0.8, 0.0,
                        0.0, 0.0, 0.1, 0.9), nrow = m, ncol = q_emiss, byrow = TRUE)
data1 <- sim_mHMM(n_t = n_t, n = n, m = m, q_emiss = q_emiss, gamma = gamma,
                  emiss_distr = emiss_distr, var_gamma = 1, var_emiss = 1)
head(data1$obs)
head(data1$states)

# including a covariate to predict (only) the transition probability matrix gamma
beta      <- rep(list(NULL), 2)
beta[[1]] <- matrix(c(0.5, 1.0,
                     -0.5, 0.5,
                      0.0, 1.0), byrow = TRUE, ncol = 2)
xx_vec      <- rep(list(NULL),2)
xx_vec[[1]] <-  c(rep(0,5), rep(1,5))
data2 <- sim_mHMM(n_t = n_t, n = n, m = m, q_emiss = q_emiss, gamma = gamma,
                  emiss_distr = emiss_distr, beta = beta, xx_vec = xx_vec,
                  var_gamma = 1, var_emiss = 1)


# simulating subject specific transition probability matrices and emission distributions only
n_t <- 0
n <- 5
m <- 3
q_emiss <- 4
gamma <- matrix(c(0.8, 0.1, 0.1,
                  0.2, 0.7, 0.1,
                  0.2, 0.2, 0.6), ncol = m, byrow = TRUE)
emiss_distr <- matrix(c(0.5, 0.5, 0.0, 0.0,
                        0.1, 0.1, 0.8, 0.0,
                        0.0, 0.0, 0.1, 0.9), nrow = m, ncol = q_emiss, byrow = TRUE)
data3 <- sim_mHMM(n_t = n_t, n = n, m = m, q_emiss = q_emiss, gamma = gamma,
                  emiss_distr = emiss_distr, var_gamma = 1, var_emiss = 1)
data3

data4 <- sim_mHMM(n_t = n_t, n = n, m = m, q_emiss = q_emiss, gamma = gamma,
                  emiss_distr = emiss_distr, var_gamma = .5, var_emiss = .5)
data4
}
\seealso{
\code{\link{mHMM}} for analyzing multilevel hidden Markov data.
}
