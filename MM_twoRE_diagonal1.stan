data {
  int<lower=1> N;
  vector[N] Y;
  int<lower=0> b;
  matrix[N,b] X;
  vector[N] D_sqrt_cov_1;  // Coviance of RE 1 is diagonal
  int<lower=1> N2;
  matrix[N2,N2] Lcov_2; // Covariance of RE 2 has rank N2 and L of LLt=Sigma is provided
}
parameters {
  vector[b] beta;
  real<lower=0> sigma;
  vector<lower=0>[2] sds;
  vector[N] z_1;
  vector[N2] z_2;
}
model {
  sigma ~ student_t(3,0,10); //half-t priors for variance component standard deviations
  sds ~ student_t(3,0,10); //half-t priors for variance component standard deviations
  z_1 ~ normal(0,1);
  z_2 ~ normal(0,1);
  Y ~ normal(X*beta + sds[1] * (D_sqrt_cov_1 .* z_1) + sds[2] * (Lcov_2 * z_2),sigma);
}
