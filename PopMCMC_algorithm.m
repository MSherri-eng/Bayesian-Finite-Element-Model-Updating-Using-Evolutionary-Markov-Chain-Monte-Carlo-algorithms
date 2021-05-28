%% Information:
% Paper Source: Evolutionary Markov Chain Monte Carlo Algorithm for Bayesian
%               model Updating.
% Code: The Population Markov Chain Monte Carlo algorithm.
% Updating five parameters of the FE model.
%--------------------------------------------------------------------------
%% The Population Markov Chain Monte Carlo algorithm.
%--------------------------------------------------------------------------
nw = length(Wc);                       %Number of updated natural frequencies
%N = ;                                 %Number of chains 
Samples = nan(nsamples,d,N);           %Generated samples matrix.
P_x = nan(T,N);                        %Matrix of chains and density
n_cr =3;                               %Defult initial value of crossover 
delta =3;                              %Defult parameters
c=0.1;                                 %Defult parameters
c_ctar=1e-12;                          %Defult parameters
p_g=0.2;                               %Defult parameters

[J,n_id]=deal(zeros(1,n_CR));           
CR=[1:n_CR]/n_CR;
pCR=ones(1,n_CR)/n_CR;                 %Crossover values and select prob. 

%Initial samples
for i=1:d
    sample(:,i)= prior(N,xmin(i),xmax(i));
end

for j=1:N
    Wc = FixedBeam_femesh(sample(j,1), sample(j,2), sample(j,3), sample(j,4), sample(j,5))';
    p_p(j,1) = Posterior(sample(j,:),Wm,Wc,Cov,Cov1,mu1); % Create initial population and compute density
    wc_x(j,1:nw)= Wc;
end

Samples(1,1:d,1:N) = reshape(sample',1,d,N); 
Wc_x(1,1:nw,1:N) = reshape(wc_x',1,nw,N);
P_x(1,1:N) = p_p'; % Store initial position of chain and density

%index of ith chain
for i = 1:N 
    
    R(i,1:N-1) = setdiff(1:N,i); 
    
end 

% Dynamic part
for t = 2:nsamples 
  
    [~,draw] = sort(rand(N-1,N)); % Randomly permute [1,...,N?1] N times
    dsample = zeros(N,d); %set N Jump vector to zero
    lambda = unifrnd(-c,c,N,1); %Draw N lambda values
    
   
    std_sample=std(sample); %calculate the std for each dimension 
    
    
    for i = 1:N
       
        D=randsample([1:delta],1,'true');  
        
        r1 = R(i,draw(1:D,i)); % Derive r1
        
        r2 = R(i,draw(D+1:2*D,i)); % Derive r2; r1 not equal r2 not equal i
        
       
        id= randsample(1:n_CR,1,'true',pCR); %select index of crossover value
       
        
        z=rand(1,d); %Draw d values from U[0,1]
        
        A = find(z<CR(id)); %drive subset of A selected dimensions.
        
        d_star = numel(A); %Determine the dimensions of samples 
        
        if d_star ==0
            [~,A]=min(z); d_star=1; %Assign at leaset one value for A 
        end
        
        gamma_d=2.38/sqrt(2*D*d_star); %Calculate the jump rate 
        
        g = randsample([gamma_d 1],1,'true',[1-p_g p_g]); %Gamma Selection
        
        dsample(i,A)=c_star*randn(1,d_star)+(1+lambda(i))* ...
            g*sum(sample(r1,A)-sample(r2,A),1);%Compute ith diff evol.
        
         
        sample_P(i,1:d)= sample(i,1:d)+dsample(i,1:d); %compute ith proposal
        
        
        
        for j = 1:d 
            overlimit_x(1,j)=sample_P(i,j)< xmax(1,j);
            underlimit_x(1,j)=sample_P(i,j)> xmin(1,j);
        end
        
        sample_P(i,1:d)=sample_P(i,1:d).*overlimit_x(1,:)+Max(1,:).*not(overlimit_x(1,:));
        sample_P(i,1:d)=sample_P(i,1:d).*underlimit_x(1,:) + Min(1,:).*not(underlimit_x(1,:));
        
        
        Wc = FixedBeam_femesh(sample_P(i,1), sample_P(i,2), sample_P(i,3), sample_P(i,4), sample_P(i,5))';
        p_px(i,1) = Posterior(sample_P(i,:),Wm,Wc,Cov,var,mu);
        wc_px(i,1:nw)=Wc;
        
        alpha = min(p_px(i,1)/p_p(i,1),1); % Compute Metropolis ratio
        
        idx = alpha > rand; % Alpha larger than U[0,1] or not?
        if idx,
            sample(i,1:d) = sample_P(i,1:d); 
            p_p(i,1) = p_px(i,1); % True: Accept proposal
            wc_x(i,1:nw)=wc_px(i,1:nw);
        else 
           dsample(i,1:d)=0;
        end 
        J(id) = J(id)+sum((dsample(i,1:d)./std_sample).^2);%Update Jump Distance Crossover
        n_id(id)=n_id(id)+1; %Corossover counter 
    
    end
    
    Samples(t,1:d,1:N) = reshape(sample',1,d,N);
    Wc_x(t,1:nw,1:N) = reshape(wc_x',1,nw,N);
    P_x(t,1:N) = p_p'; % Add current position and density to chain
    
    if t<T/10
        pCR=J./n_id;
        pCR=pCR/sum(pCR); %Update the selection prob.
    end 
 
   %---------------------------------------------
  % Calculate the updated values and the TAE.
    samples_1 = Samples(:,:,1);
    Ehcv = (1/t)*sum(samples_1,1);
    Wccv = FixedBeam_femesh(Ehcv(1),Ehcv(2), Ehcv(3), Ehcv(4), Ehcv(5))';

    Percentagecv = zeros(1,nw);

    for j = 1:nw
        Percentagecv(1,j) = 100*(abs(Wccv(j) - Wm(j))/Wm(j));
    end
    
    Errorcv = [Errorcv, (sum(Percentagecv))/nw];
    
end
