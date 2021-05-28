%% Information:
% Paper Source: Evolutionary Markov Chain Monte Carlo Algorithm for Bayesian
%               model Updating.
% Code: Bayesian and updating parameters modelling used for updating algorithm
% Updating five parameters of the FE model.
%--------------------------------------------------------------------------
%% Bayesian and updating parameters modelling 
%--------------------------------------------------------------------------
% Posterior function:
%Wm:    measured natural frequencies. 
%Wc:    calculated natural frequencies.
%sup:   sampled updating parameter value.
%mup:   mean updating parameter value.
%var:   varince vector.
%Cov:   covaiance matrxi. 
Posterior = @(Wm,Wc,Cov,var,mu,sup) (Wc - Wm)*(inv(var*var))*(Wc - Wm)'...
    +(sup-mup)*(inv(Cov*Cov))*(sup-mup)';
%--------------------------------------------------------------------------
% Prior function:
%x_par:     updating parameters.
%sig:       step control for single-chain method.
%d:         dimensional vector of the updating parameters.
%N = ;         number of chains for multichain algorithms. 
Pior = @(x_par,mu,sig,d,N) (1/(2*pi)^(d/2))*(1/det(sig*sig)^0.5)...
    *exp(-0.5*(x_par-mu)*(inv(sig*sig))*(x_par-mu)');
%--------------------------------------------------------------------------
%The initial 
%nsamples = ;                          %number of required samples.
%uE_m = ;                              %updating element mean value.
%uE_c = ;                              %updating element calculated value.
Par1 = uE_c;                           %updating parameter number 1.
Par2 = uE_c;                           %updating parameter number 2.
Par3 = uE_c;                           %updating parameter number 3.
Par4 = uE_c;                           %updating parameter number 4.
Par5 = uE_c;                           %updating parameter number 5.
x_par = [Par1, Par2, Par3, Par4, Par5];
d = length(x_par);
%Wm = [] ;
%--------------------------------------------------------------------------
%Call FEA simulator to calculate the natural frequencies
Wc = FEA_simulator(Par1, Par2, Par3, Par4, Par5)';
W_init = Wc;                           %Initial calculated natural frequencies.
%--------------------------------------------------------------------------
%% Call the updating algorithm 
Errorcv = [];                          %Matrix of total evaluation error.                               
samples = zeros(nsamples, d);          %Matrix of returned samples.
% xmax = [];                           %Maximum bound vector of d_parameters
% xmin = [];                           %Minimum bound vector of d_parameters
blength = xmax-xmin;                   %Bound lenght
%--------------------------------------------------------------------------
%% Results and Figures 
Par_f = (1/nsamples)*sum(samples,1);        %mean values of the updated parameters
Wc_f = FEA_simulator(Par_f(1),Par_f(2), Par_f(3),Par_f(4),Par_f(5))'; %Final Wc.
StandardD = std(samples,0,1);               %The Standard Deviation of the samples
corr = cov(samples)./(StandardD'*StandardD);%The correlation between the updated parameters
figue(); plot(Errorcv);                     %Plot error curve 
figure(); bar3(corr);                       %Plot the corr.
figure(); boxplot(samples);                 %box-plot of the generated samples
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------


