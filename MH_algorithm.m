%% Information:
% Paper Source: Evolutionary Markov Chain Monte Carlo Algorithm for Bayesian
%               model Updating.
% Code: The Metropolis-Hastings algorithm.
% Updating five parameters of the FE model.
%--------------------------------------------------------------------------
%% The Metropolis-Hastings algorithm.
%--------------------------------------------------------------------------
logpOld = Posterior(x_par,Wm,Wc,Cov,var,mu);       %Posterior of Initial Paramters
covsamples = [];                                   %Generated samples
%The dynamic part
for t=1:nsamples
   
    xprime = x_par + randn(1,d)*sig;
    
    for i = 1:d 
      overlimit_x(1,i)=xprime(1,i)<=xmax(1,i);
      underlimit_x(1,i)=xprime(1,i)>=xmin(1,i);
    end
  
   RV_interval = rand(1,d).*blength + xmin;
   xprime=xprime(1,:).*overlimit_x(1,:)+xmax(1,:).*not(overlimit_x(1,:));
   xprime=xprime(1,:).*underlimit_x(1,:) + xmin(1,:).*not(underlimit_x(1,:));
   Wc = FEA_simulator(xprime(1,1), xprime(1,2),xprime(1,3), xprime(1,4), xprime(1,5))';
   logpNew = Posterior(xprime,Wm,Wc,Cov,var,mu);     
   alpha = exp(logpOld-logpNew);     
   qnumer = proposalProb(x_par,xprime,sig,d);   %proposalProb; % q(x|x')         
   qdenom = proposalProb(xprime,x_par,sig,d);   %proposalProb; % q(x'|x)         
   alpha = alpha * (qnumer/qdenom);
   r = min(1, alpha);
   u = rand(1,1);
    if u < r         
       x_par = xprime;         
       naccept = naccept + 1;         
       logpOld = logpNew; 
       Prob1 = r;
    end
   Prob = [Prob; Prob1];
   samples(t,:) = x_par;
   Wc = FEA_simulator(x_par(1,1), x_par(1,2),x(1,3), x_par(1,4),x_par(1,5))';
   WWc(t,:) = Wc;
   covsamples = [covsamples; x];
   Ehcv = (1/t)*sum(covsamples,1);      %Current updated parameter values.                           
   Wccv = FEA_simulator(Ehcv(1),Ehcv(2), Ehcv(3),Ehcv(4), Ehcv(5))'; %Wc current value.
   Percentagecv = zeros(1,length(Wc));  %Percentage error. 
    for j = 1:nparams
        Percentagecv(1,j) = 100*(abs(Wccv(j) - Wm(j))/Wm(j));
    end
   Errorcv = [Errorcv, (sum(Percentagecv))/length(Wc)]
 end
