%% Information:
% Paper Source: Evolutionary Markov Chain Monte Carlo Algorithm for Bayesian
%               model Updating.
% Code: The Differential Evolution Markov Chain Monte Carlo algorithm.
% Updating five parameters of the FE model.
%--------------------------------------------------------------------------
%% The Differential Evolution Markov Chain Monte Carlo algorithm.
%--------------------------------------------------------------------------
nw = length(Wc);                       %Number of updated natural frequencies
gamma_RWM = 2.38/sqrt(2*d);            %dEFULT jump rate
%N = ;                                 %Number of chains 
Samples = nan(nsamples,d,N);           %Generated samples matrix.
Samples(1,1:d,1:N) = reshape(sample',1,d,N); 

%Initial samples
for i=1:d
    sample(:,i)= prior(N,xmin(i),xmax(i));
end

%Initial calculated natural frequencies 
for j=1:N
    Wc = FEA_simulator(sample(j,1), sample(j,2), sample(j,3), sample(j,4), sample(j,5))';
    I_p(j,1) = Posterior(sample(j,:),Wm,Wc,Cov,var,mu); %  initial population 
    wc_x(j,1:nw)= Wc;
end

%the index of chains for DE
for i = 1:N 
    R(i,1:N-1) = setdiff(1:N,i); 
end 

% Dynamic part
for t = 2:nsamples 
    [~,draw] = sort(rand(N-1,N)); 
        g = randsample([gamma_RWM 1],1,true,[0.9 0.1]); 
    for i = 1:N 
        r1 = R(i,draw(1,i)); 
        r2 = R(i,draw(2,i)); 
        sample_P(i,1:d) = sample(i,1:d) + g*(sample(r1,1:d)-sample(r2,1:d))...
            + (1e-6)*randn(1,d);
        
        for j = 1:d 
            overlimit_x(1,j)=sample_P(i,j)< xmax(1,j);
            underlimit_x(1,j)=sample_P(i,j)> xmin(1,j);
        end
        
        sample_P(i,1:d)=sample_P(i,1:d).*overlimit_x(1,:)+Max(1,:).*not(overlimit_x(1,:));
        sample_P(i,1:d)=sample_P(i,1:d).*underlimit_x(1,:) + Min(1,:).*not(underlimit_x(1,:));
        
        
        Wc = FEA_simulator(sample_P(i,1), sample_P(i,2), sample_P(i,3), sample_P(i,4), sample_P(i,5))';
        p_px(i,1) = Posterior(sample_P(i,:),Wm,Wc,Cov,var,mu);
        wc_px(i,1:nw)=Wc;
        
        alpha = min(p_px(i,1)/p_p(i,1),1); 
        
        idx = alpha > rand; 
        if idx,
            sample(i,1:d) = sample_P(i,1:d); 
            p_p(i,1) = p_px(i,1); 
            wc_x(i,1:nw)=wc_px(i,1:nw);
        end
    end
    
    Samples(t,1:d,1:N) = reshape(sample',1,d,N);
    Wc_x(t,1:nw,1:N) = reshape(wc_x',1,nw,N);
    P_x(t,1:N) = p_p'; 
    
    samples_1 = Samples(:,:,1);
    samples_1(isnan(samples_1)) = 0.0;
    
    Ehcv = (1/t)*sum(samples_1,1);
    Wccv = FEA_simulator(Ehcv(1),Ehcv(2), Ehcv(3), Ehcv(4), Ehcv(5))';

    Percentagecv = zeros(1,nw);

    for j = 1:nw
        Percentagecv(1,j) = 100*(abs(Wccv(j) - Wm(j))/Wm(j));
    end
    
    Errorcv = [Errorcv, (sum(Percentagecv))/nw];
    
end