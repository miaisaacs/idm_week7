## SIR MODEL
## SCRIPT TO TEST THE SIR MODEL run deterministically vs stochastically

## if you haven't had it installed yet, type in the following first:
## install.packages('deSolve')
library(deSolve)

################################################################################
## SIMULAITON USING THE SIR MODEL
# General Routine: 
# Step 1 - Code the ODE/model (or function here) [note the SIR model is described by a set of ODEs]
# Step 2 - Specify the initial conditions/parameters
# Step 3 - Call the “ode” function to solve the ODE/model and generate the simulation
# Step 4 - Check/analyze model outputs
################################################################################
# Step 1: code the SIR model (or function here)
SIR=function(t,state,parameters){
  with(as.list(c(state,parameters)),{
    # rate of change
    dS = -beta*S*I/N;
    dI= beta*S*I/N - gamma*I;
    
    # Optional: also compute the cumulative number of infection
    dcumI = beta*S*I/N
    
    # return the rate of change
    list(c(dS, dI, dcumI))
  }) # end with(as.list...)
}

SIRstoch=function(t,state,parameters){
  with(as.list(c(state,parameters)),{
    # get parameters (fixed) for each time step
    beta_t = beta[floor(t)+1] 
    gamma_t = gamma[floor(t)+1] 
    
    # rate of change
    dS = -beta_t *S*I/N;
    dI= beta_t *S*I/N - gamma_t *I;
    # Optional: also compute the cumulative number of infection
    dcumI = beta_t *S*I/N
    
    # return the rate of change
    list(c(dS, dI, dcumI))
  }) # end with(as.list...)
}
# Step 2: specify initial conditions/parameters
N=1e5; # population size
I0=10; # initial No. of Infectious people
S0=N-I0;  # initial No. of Susceptible people
state=c(S=S0, I=I0, cumI = I0);  # store the inital conditions I & S together 
beta = .5; gamma = .3
parameters=c(beta=beta,gamma=gamma);  # store the model parameters, unit: per day

times=seq(1,100,by=1);  # set simulation time steps: 1 to 100 days here

# Step 3: call the ode function to generate the simulation
sim=ode(y=state,times=times,func=SIR,parms=parameters);

# n stochastic runs
n_runs = 100
S_stoch = I_stoch = cumI_stoch = matrix(0, nrow = length(times), ncol=n_runs)
for(i in 1:n_runs){
  # run it with the stochastic SIR model
  # add noise to key parameters, e.g., here small Gaussian noise mean = 0, SD = .05
  parameters=list(beta= beta * (1 + rnorm(length(times)+1, 0, 0.05)), 
                  gamma= gamma * (1 + rnorm(length(times)+1, 0, 0.05)));  
  
  sim_stoch_t = ode(y=state,times=times,func=SIRstoch,parms=parameters);
  # save the results
  S_stoch[,i] = sim_stoch_t[,'S']
  I_stoch[,i] = sim_stoch_t[,'I']
  cumI_stoch[,i] = sim_stoch_t[,'cumI']
  print(i)
}

# Step 4: check model outputs
# compare deterministic SIR with n stochastic runs
matplot(S_stoch, col='grey', type='l', xlab='time', ylab='S')
lines(sim[,'S'])

matplot(I_stoch, col='grey', type='l', xlab='time', ylab='I')
lines(sim[,'I'])

matplot(cumI_stoch, col='grey', type='l', xlab='time', ylab='cumI')
lines(sim[,'cumI'])


