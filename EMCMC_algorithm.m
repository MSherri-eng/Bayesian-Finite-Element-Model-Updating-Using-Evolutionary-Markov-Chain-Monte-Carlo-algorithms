%% Information:
% Paper Source: Evolutionary Markov Chain Monte Carlo Algorithm for Bayesian
%               model Updating.
% Code: The Evolutionry Markov Chain Monte Carlo algorithm.
% Updating five parameters of the FE model.
%--------------------------------------------------------------------------
%% The Evolutionry Markov Chain Monte Carlo algorithm.
nw = length(Wc);   % Number of updated natural frequencies
%MaxSubIt = ;      % Maximum Number of Sub-iterations
%T0= ;             % Initial Temp.
%alpha= ;          % Parallel temp. reduction rate
%nPop = ;          % Population Size
%nMove =  ;        % Number of Moves

% Create Empty Structure for Individuals
empty_individual.Position=[];
empty_individual.Cost=[];

% Create Population Array
pop=repmat(empty_individual,nPop,1);

% Initialize Best Solution
BestSol.Cost=inf;


% Initialize Population
for i=1:nPop
    
    % Initialize Position
     mx=xmax(1);       %Max bound
     mn=xmin(1);       %Max bound

     Max = [mx,mx,mx,mx,mx];
     Min = [mn,mn,mn,mn,mn];

     n=model.n;
    
     %Initial random positions
     sol=Min+rand(1,n).*(Max-Min); 
     pop(i).Position=sol;
    
    % Evaluation
    Wc= FEA_simulator(pop(i).Position(:,1), pop(i).Position(:,2), pop(i).Position(:,3),pop(i).Position(:,4), pop(i).Position(:,5))';
    pop(i).Cost=Posterior(pop(i).Position,Wm,Wc,Cov,var,mu);
    
    
    % Update Best Solution
    if pop(i).Cost<=BestSol.Cost
        BestSol=pop(i);
    end
    
end

% Array to Hold Best Cost Values
BestCost=zeros(nsamples,1);

% Intialize Temp.
T=T0;

%Dynamic Part
for it=1:nsamples
    
    for subit=1:MaxSubIt
        
        % Create and Evaluate New Solutions
        newpop=repmat(empty_individual,nPop,nMove);
        for i=1:nPop
            for j=1:nMove
                
              % Create Neighbor
                   pSwap=0.2;          %Defult set
                   pReversion=0.5;     %Defult set
                   pInsertion=1-pSwap-pReversion;

                   p=[pSwap pReversion pInsertion];

               %selection 
                   r=rand;
                   c=cumsum(p);
                   i=find(r<=c,1,'first'); 
                        
                 switch 
                   case 1
                   % Swap
                   n=numel(Position1);
                   I=randsample(n,2);
                   i1=I(1);
                   i2=I(2);
                   Position2([i1 i2])=Position1([i2 i1]);
          
                   case 2
                   % Reversion
                   n=numel(Position1);
                   I=randsample(n,2);
                   i1=min(I);
                   i2=max(I);
                   Position2=Position1;
                   Position2(i1:i2)=Position1(i2:-1:i1);

                   case 3
                   % Insertion
                    n=numel(Psition1);
                    I=randsample(n,2);
                    i1=I(1);
                    i2=I(2);
                if i1<i2
                Position2=Positio1([1:i1-1 i1+1:i2 i1 i2+1:end]);
                else
                Psition2=Position1([1:i2 i1 i2+1:i1-1 i1+1:end]);
                end
                                
      Position=Position2;
      newpop(i,j).Position= pop(i).Position);
                
     for d = 1:nparams 
       overlimit_x(1,d)=newpop(i,j).Position(d)<=xmax(1,d);
       underlimit_x(1,d)=newpop(i,j).Position(d)>=xmin(1,d);
     end

    RV_interval = rand(1,nparams).*blength + Min;
    newpop(i,j).Position=newpop(i,j).Position(1,:).*overlimit_x(1,:)+RV_interval.*not(overlimit_x(1,:));
    newpop(i,j).Position=newpop(i,j).Position(1,:).*underlimit_x(1,:) + RV_interval.*not(underlimit_x(1,:));
                
   % Evaluation
   Wc= FEA_simulator(newpop(i,j).Position(:,1), newpop(i,j).Position(:,2), newpop(i,j).Position(:,3),newpop(i,j).Position(:,4), newpop(i,j).Position(:,5))';
   newpop(i,j).Cost=Posterior(newpop(i,j).Position,Wm,Wc,Cov,Cov1,mu1);
        
   end
  end
newpop=newpop(:);
        
% Sort Neighbors
[~, SortOrder]=sort([newpop.Cost]);
newpop=newpop(SortOrder);
        
   for i=1:nPop
            
      if newpop(i).Cost<=pop(i).Cost
       pop(i)=newpop(i);
                
       else
         DELTA=(newpop(i).Cost-pop(i).Cost)/pop(i).Cost;
          P=exp(-DELTA/T);
          if rand<=P
           pop(i)=newpop(i);
                
         end
        end
            
            % Update Best Solution Ever Found
            if pop(i).Cost<=BestSol.Cost
                BestSol=pop(i);
                
            else
                %samples(end+1,:)= BestSol.Position+rand(1,nparams);
            end
        
        end

    end
   
   
    % Store Best Cost Ever Found
    BestCost(it)=BestSol.Cost;
    samples(end+1,:)=BestSol.Position;
    
    
    % Display Iteration Information
    disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(BestCost(it))]);
    
    % Update Temp.
    T=alpha*T;
    

 
    Ehcv = (1/it)*sum(samples,1);
    
    Wccv = FEA_simulator(Ehcv(1),Ehcv(2), Ehcv(3), Ehcv(4), Ehcv(5))';

    Percentagecv = zeros(1,nw);

    for p = 1:nw
        Percentagecv(1,p) = 100*(abs(Wccv(p) - Wm(p))/Wm(p));
    end
    
    Errorcv = [Errorcv, (sum(Percentagecv))/nw]
   
    
end