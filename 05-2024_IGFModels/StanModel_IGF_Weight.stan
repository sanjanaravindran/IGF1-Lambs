data {
  int<lower=0> N; //number of data points
  //vector[N] StorageTime_sc; // storage time predictor
  vector[N] Weight; //bodyweight response
  vector[N] IGF1_sc; //igf predictor
  vector[N] PopSize_sc; //pop size predictor
  vector[N] SexF; // sex predictor
  vector[N] Twin; // twin predictor
  int<lower=1> num_years; //number of years
  int<lower=1, upper=num_years> BirthYear[N]; //year id
  int<lower=1> num_mums; //number of mums
  int<lower=1, upper=num_mums> MumID[N]; //mum id
  int<lower=1> num_days; //number of days
  int<lower=1, upper=num_days> ELISARunDate[N]; //day id
  int<lower=1> num_plates; //number of plates
  int<lower=1, upper=num_plates> PlateNumber[N]; //plate id
}

parameters {
  //real beta_ST; // slope for storage time
  real alpha; //intercept for weight~IGF + ... model
  real beta_IGF; //slope for IGF in weight~IGF + ... model
  real beta_SexF; // slope for sex in weight~IGF + ... model
  real beta_Twin; // slope for twin in weight~IGF + ... model
  real beta_PopSize; // slope for pop size in weight~IGF + ... model
  vector[num_years] t; //year intercepts
  vector[num_mums] u; //mum intercepts
  vector[num_days] v; //day intercepts
  vector[num_plates] w; //plate intercepts
  real<lower=0> sigma_e; //error sd
  real<lower=0> sigma_t; //year sd
  real<lower=0> sigma_u; //year sd
  real<lower=0> sigma_v; //day sd 
  real<lower=0> sigma_w; //plate sd 
  real<lower=0> sigma_e1;
  
  
}

transformed parameters {
  vector[N] IGF_true = IGF1_sc - v[ELISARunDate] - w[PlateNumber]; 
}

model {
  sigma_e ~ cauchy(0,2.5); //prior scale
  sigma_e1 ~ cauchy(0,2.5);  //prior scale
  
  
  //beta_IGF ~ std_normal(); //prior location
  //beta_SexF ~ std_normal(); //prior location
  //beta_Twin ~ std_normal(); //prior location
  //beta_PopSize ~ std_normal(); //prior location 
  
  //NCP stuff
  t ~ normal(0,1);
  u ~ normal(0,1);
  v ~ normal(0,sigma_v);
  w ~ normal(0,sigma_w);
  sigma_t ~ cauchy(0,2.5);
  sigma_u ~ cauchy(0,2.5);
  sigma_v ~ cauchy(0,2.5);
  sigma_w ~ cauchy(0,2.5);
  
  //likelihood
  IGF1_sc ~ normal(v[ELISARunDate] + w[PlateNumber], sigma_e1);
  
  // likelihood
  Weight ~ normal(alpha + beta_IGF * IGF_true + beta_SexF * SexF + 
                    beta_Twin * Twin + beta_PopSize * PopSize_sc +
                    sigma_t*t[BirthYear] + sigma_u*u[MumID], sigma_e);
}

generated quantities {
  vector[N] Weight_rep;
  for (i in 1:N)
    Weight_rep[i] = normal_rng(alpha  + beta_IGF * IGF_true[i] + beta_SexF * SexF[i] + 
                                 beta_Twin * Twin[i] + beta_PopSize * PopSize_sc[i] + 
                                 sigma_t*t[BirthYear[i]] + sigma_u*u[MumID[i]], sigma_e);
}
