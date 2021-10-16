%-------------------------------------------------------------------------%
%                      AIYAGARI (1994) MODEL                              %
%-------------------------------------------------------------------------%
% Reference: Slides of Advanced Macroeconomics II, University of Konstanz.
% Course material available at:
% https://github.com/aledinola/Advanced-Macroeconomics-II

% Household maximization problem:
% V(z) = max_{a_hat'} {u(z-a_hat')+beta*EV(z')}
% s.t.
% c+a_hat'=z, z'=(1+r)a_hat'+w*s-r*phi, a_hat'>=0
% where z is cash on hand, a_hat=a+phi. phi can be either zero (ad hoc)or
% ws(1)/r (natural borrowing limit). See slides for more details.
% Marginal factor pricing:
% r = F_1(K,L)-delta, wage=F_2(K,L)
% Market clearing:
% K = integral of a
% -------------------------------------------------------------------------
% Original source: Leo Kaas GAUSS code. Translated in MATLAB by Alessandro
% Di Nola.
% -------------------------------------------------------------------------

clear
clc
close all

disp('SOLVING THE AIYAGARI MODEL')

% Set some flags
par.do_plots = 0;
par.verbose  = 1;

% Computational parameters
demtol  = .001;   % tolerance for equilibrium condition 
par.eps = .00001; % tolerance for value function iteration
par.concavity = 0; % assume value function is concave

% Initialize model parameters

par.alpha = .3;  % capital coefficient in production function
par.beta  = .95; % discount factor
par.gam   = 3;   % CRRA parameter
par.del   = .1;  % capital depreciation rate

% Grid for resources z

zdim = 1000; % grid size 
zmin = .01;  % Set zmin = very small number>0 to avoid numerical problems
zmax = 50;

% Recall to set a(1)=0, a(N)=z(N-1)
z_grid = linspace(zmin,zmax,zdim)'; % vector of resources (state variable)
a_grid = [0;z_grid(1:end-1)];       % vector of assets (a-hat, choice variable)

% Labor supply vector

sdim = 2;    % size of idiosyncratic productivity state space
smin = .2;   % range for productivities
smax = 1.8;

s_grid  = linspace(smin,smax,sdim)'; % productivity grid */
pi_s = (1/sdim)*ones(sdim,1);     % probabilities */

% IMPORTANT: Note that E(s)=1 This is also the amount of aggregate labor.

%% Start general equilibrium loop
rr0 = [0.005,0.03]; % Initial condition for the interest rate

[rr,FVAL,EXITFLAG] = fzero(@(x) excess_demand(x,a_grid,z_grid,s_grid,pi_s,par),rr0);
    
[ED,Ea,sol,sim,phi] = excess_demand(rr,a_grid,z_grid,s_grid,pi_s,par);

save results

%% Save results in txt file

if phi==0
    fileID = fopen('results_nobor.txt','w');
    fprintf(fileID,'Borrowing limit   = %12.8f \n',-phi);
    fprintf(fileID,'Interest rate     = %12.8f \n',rr);
    fprintf(fileID,'Aggregate capital = %12.8f \n',Ea);
    fprintf(fileID,'Excess demand     = %12.8f \n',ED);
    fclose(fileID);
else
    fileID = fopen('results_naturalbor.txt','w');
    fprintf(fileID,'Borrowing limit   = %12.8f \n',-phi);
    fprintf(fileID,'Interest rate     = %12.8f \n',rr);
    fprintf(fileID,'Aggregate capital = %12.8f \n',Ea);
    fprintf(fileID,'Excess demand     = %12.8f \n',ED);
    fclose(fileID);
end


%% Plot some results
histogram(sim.asset_sim,'Normalization','probability')
xlabel('Assets, a')
ylabel('Density')
if phi==0
    title('Asset distribution with \Phi=0')
    print('asset_dist_nobor','-dpng')
else
    title('Asset distribution with \Phi=s_1*w/r')
    print('asset_dist_naturalbor','-dpng')
end


figure
subplot(3,1,1)
plot(sim.s_sim(end-300:end),'linewidth',2)
xlim([0,300])
title('Productivity shock, s')
subplot(3,1,2)
plot(sim.asset_sim(end-300:end),'linewidth',2)
xlim([0,300])
title('Assets, a')
subplot(3,1,3)
plot(sim.c_sim(end-300:end),'linewidth',2)
xlim([0,300])
title('Consumption, c')
xlabel('Time (simulation periods)')
if phi==0
    print('simulation_nobor','-dpng')
else
    print('simulation_naturalbor','-dpng')
end





