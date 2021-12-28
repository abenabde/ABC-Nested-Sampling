function theta = abc_sdof(data,tol,accuracy)
%::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
% This file demonstrates the use of the Approximate 
% Bayesian Computation based on Nested Sampling (ABC-NS)
% algorithm to deal with parameter estimation
% 
% The training data is simulated from a SDOF system 
%                  mx'' + cx' + kx = f 
% m : Mass of the system (supposed to be known) 
% c : Damping coefficient
% k : Stiffness coefficient
% f : Harmonic excitation 
%
%::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
global m c k
ts = cputime;
Av = 100/(length(data)*std(data)^2); % normalisation constant
n=length(data); % number of observations in the training data
beta_0=0.6; % proportion of the remaining ('alive') particles 
alpha_0=0.3; % proportion of the dropped particles
f_0=1.1; % enlargement factor
prior_c = [10 30]; % uniform prior on damping parameter
prior_k = [50 200]; % uniform prior on stiffness parameter 
n_sim=100; % number of ABC-NS iterations
n_param=2; % number of parameters to be inferred
threshold_vec(1,1) = tol; % set the initial value of the tolerance threshold
pop_size = 1000; % population size
% create matrix to store the particles for all the populations
population = ones(pop_size,n_sim*n_param); 
nmse_val = ones(pop_size,n_sim*1);
acceptance = zeros(n_sim,1); % store the acceptance rate
weights = ones(pop_size,n_sim); % weights matrix
accuracy; %  The convergence precision
iter = 0;
sim_0 = 0;
sim=1; % population number
%-----------------------------------------------------------
%  STEP 1 : First loop/ Get the first population
%-----------------------------------------------------------
fprintf('~~~~~~~~~~~~~~~~~~~~ABC-NS algorithm: Simulation has started ~~~~~~~~~~~~~~~~~~~~\n ')
while iter < pop_size
       sim_0 = sim_0 + 1;
       % Generate a particle from the priors
          theta_2s = unifrnd([prior_c(1,1) prior_k(1,1)],[prior_c(1,2) prior_k(1,2)]);        
          % Simulate from the model 
            [t,x0] = model(theta_2s); 
           %  Compute the Normalised Mean Square Error (NMSE) (used to measure the discrepancy)
           % between the observed and the simulated data
            NMSE = Av*sum((data - x0(:,1)).^2);               
           if NMSE <= threshold_vec(1,1)
              pop(iter+1,:) = theta_2s;
              nmse_v(iter+1,:) =  NMSE;
              iter = iter +1;
           end           
end  
population(:,1:n_param) = pop; % matrix to store the particles 
acceptance(1,1) = pop_size/sim_0; % compute the acceptance rate
[svc,ind] = sort(nmse_v,'descend');
threshold_vec(1,2) = svc(alpha_0*pop_size);  % Define the next tolerance threshold value 
% Assign a weight to each particle
for i=1:pop_size      
       if nmse_v(i,1)>=threshold_vec(1,2) 
           weight(i,:)=0;
       else 
           weight(i,:)=(1/threshold_vec(1,1))*(1 - (nmse_v(i,1)./threshold_vec(1,1)).^2);
       end
end
[seq,idx] = datasample([pop],pop_size*beta_0,'Weights',[weight'],'Replace',false);
active_pop = seq;
nmse_pop = nmse_v(idx);
% create matrix to store the remaining 'alive' particles
active_samples = ones(pop_size*beta_0,n_param*n_sim); 
% create a matrix to store the NMSE values associated to the remaining particles
active_obj = ones(pop_size*beta_0,n_sim); 
active_samples(:,1:n_param) = active_pop; % store the 'alive' particles for the first population
active_obj(:,1) = nmse_pop;  % store the associated NMSE values
%--------------------------------------------------------------------------
% Determine the optimum ellipse 
mu = mean(active_pop); % mass center of the remaining particles
B = cov(active_pop); % covariance matrix based on the remaining particles
D=n_param; % number of parameters to be identified
% calculate volume of bounding ellipsoid
const = pi^(D/2)/gamma(D/2 + 1);
VS = const*sqrt(det(B));
[Bs, mus, VEs, ns] = calc_ellipsoid(active_pop, VS);
Bsn = Bs*f_0;
fprintf('\nABC-NS population #%d --- Tolerance threshold %4.3f',sim,threshold_vec(1,sim))
% %========================================================================
%                  STEP 2 : Second loop/ Get next populations
% %========================================================================
for sim=2:n_sim
iter =0;
weight=[];
active_pop=[];
nmse_pop=[];svc=[];
next_pop = [];comp_pop =[];nmse_values=[];
idx=[];mu=[];B=[];Vs=[];Bs=[];VEs=[];sin0=0;
while iter<pop_size*0.4
             sm=0;
             while sm<1
                 % Sample a new particle inside the ellipse
             theta_new = draw_from_ellipsoid(Bsn, mus, 1);       
               %----------------Uniform Kernel ----------------------------
              if (theta_new(1,1)>=min(seq(:,1)) && theta_new(1,1)<=max(seq(:,1)) ...
                      && theta_new(1,2)>=min(seq(:,2)) &&  theta_new(1,2)<=max(seq(:,2)))  
                   sm = sm + 1;
                   break
               end                              
             end
            % Simulate data 
            [t,x0] = model(theta_new);
           % Compute the NMSE 
            NMSE = Av*sum((data - x0(:,1)).^2);  
               if  NMSE <= threshold_vec(1,sim)
                   comp_pop(iter+1,:) = theta_new;
                   nmse_values(iter+1,:) =  NMSE;
                   iter = iter +1 ;
               end     
            sin0 = sin0 + 1 ;             
end
           % Replenish the population by replacing the dropped particles
               next_pop = [active_samples(:,n_param*sim-(2*n_param-1):n_param*sim-n_param);comp_pop];
               nmse_store =[active_obj(:,sim-1);nmse_values];
               [svc,ind] = sort(nmse_store,'descend');         
               % Compute the acceptance rate 
               acceptance(sim,1) = (pop_size*0.4)/sin0;
               % Compute the next tolerance threshold value
                threshold_vec(1,sim+1) = svc(alpha_0*pop_size) ;
      %----------------------------------                   
% Assign a weight for the new particles
for vv=1:1000    
       if nmse_store(vv,1)>=threshold_vec(1,sim)
           weight(vv,:)=0;
       else 
           weight(vv,:)=(1/threshold_vec(1,sim))*(1 - (nmse_store(vv,1)./threshold_vec(1,sim)).^2);
       end
end
% Allows to select the remaining particles in a probabilistic way
               [seq,idx] = datasample([next_pop],pop_size*beta_0,'Weights',[weight],'Replace',false);
               active_pop = seq;
               nmse_pop = nmse_store(idx);
   %-------------------------------
   % store active samples and the associated NMSE  
               active_samples(:,n_param*sim-(n_param-1):n_param*sim) = active_pop;
               active_obj(:,sim) = nmse_pop;   
         % Update the mass center and the covaraince matrix
            mu = mean(active_pop); % mass center
            B = cov(active_pop); % covariance matrix
            % calculate volume of bounding ellipsoid
            const = pi^(D/2)/gamma(D/2 + 1);
            VS = const*sqrt(det(B)); 
           Bsn = [];mus=[];
          % Get Bs and mus
            [Bs, mus, VEs, ns] = calc_ellipsoid(active_pop, VS);
          % Enlarge the ellipse 
            Bsn = Bs*f_0;
      %---------------------------------------------------- 
      % store the particles and the associated NMSE values
                population(:,n_param*sim-(n_param-1):n_param*sim) = [active_samples(:,n_param*sim-(2*n_param-1):n_param*sim-n_param) ;comp_pop];
                nmse_val(:,sim) = [active_obj(:,sim);nmse_values];
%          Stopping criterion 
                if abs(threshold_vec(1,sim,1) - threshold_vec(1,sim + 1))<accuracy
                    break
                end    
           % Notice we could decrease further this value (accuracy),
           % of course the CPU time will increase

  fprintf('\nABC-NS population #%d --- Tolerance threshold %4.3f',sim,threshold_vec(1,sim))
end 
fprintf('\n~~~~~~~~~~~~~~~~~~~~ABC-NS algorithm has converged ~~~~~~~~~~~~~~~~~~~~\n ')
theta=mean(next_pop);
figure(1)
subplot(2,2,1)
hist(next_pop(:,1))
xlabel('c');ylabel('Frequency')
subplot(2,2,2)
hist(next_pop(:,2))
xlabel('k');ylabel('Frequency')
subplot(2,2,3:4)
plot(next_pop(:,1),next_pop(:,2),'k.')
set(gcf,'Color','w')
axis([10 30 50 200])
xlabel('c');ylabel('k')
cpu = (cputime - ts)/60;
fprintf('\n~~~~~~~~~~~~~~~~~~~~ Estimated CPU : %2.4f minutes', cpu) 
% Activate to save and postprocess the results
% fname = sprintf('ABCNS_SDOF.mat', sim);
% save(fname) 