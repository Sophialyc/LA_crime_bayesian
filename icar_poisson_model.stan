functions {
    real icar_normal_lpdf(vector phi, int N, array[] int node1, array[] int node2) {
        return -0.5 * dot_self(phi[node1] - phi[node2]);
    }
}
data {
    int<lower=0> N;
    int<lower=0> N_edges;
    array[N_edges] int<lower=1, upper=N> node1;
    array[N_edges] int<lower=1, upper=N> node2;
    array[N] int<lower=0> Y;                                        // dependent variable i.e., number of road accidents
    vector<lower=0>[N] X;                                           // independent variable i.e., deprivation score
    vector<lower=0>[N] E;                                           // estimated number of expected cases of road accidents
}
transformed data {
    vector[N] log_offset = log(E);                                  // use the expected cases as an offset and add to the regression model
}
parameters {
    real alpha;                                                     // define the intercept (overall risk in population)
    real beta;                                                      // define the coefficient for the deprivation score variable
    real<lower=0> sigma;                                            // define the overall standard deviation producted with spatial effect smoothing term phi
    vector[N] phi;                                                  // spatial effect smoothing term or spatial ICAR component of the model 
}
model {
    phi ~ icar_normal(N, node1, node2);                             // prior for the spatial random effects
    Y ~ poisson_log(log_offset + alpha + beta*X + phi*sigma);       // likelihood function i.e., spatial ICAR model using Possion distribution
    alpha ~ normal(0.0, 1.0);                                       // prior for intercept   (weak/uninformative prior)
    beta ~ normal(0.0, 1.0);                                        // prior for coefficient (weak/uninformative prior)
    sigma ~ normal(0.0, 1.0);                                       // prior for SD          (weak/uninformative prior)
    sum(phi) ~ normal(0, 0.001*N);
}
generated quantities {
    vector[N] eta = alpha + beta*X + phi*sigma;                     // do eta equals alpha + beta*X + phi*sigma to get the relative risk for areas 
    vector[N] mu = exp(eta);                                        // the exponentiate eta to mu areas-specific relative risk ratios (RRs)
}

