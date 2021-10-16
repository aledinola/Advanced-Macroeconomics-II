function [ED,Ea,sol,sim,phi] = excess_demand(rr,a_grid,z_grid,s_grid,pi_s,par)

% Unpack structure
del   = par.del;
alpha = par.alpha;
beta  = par.beta;
gam   = par.gam;
eps   = par.eps;
verbose = par.verbose;
do_plots = par.do_plots;
concavity = par.concavity;

zdim = size(z_grid,1);
sdim = size(s_grid,1);

Kd = capital_demand(rr,del,alpha);
ww = wage(Kd,alpha);

% Debt limit (two options)
phi = 0;               % zero borrowing constraint
% if rr>0
%     phi = ww*s_grid(1)/rr; % natural debt limit
% else
%     phi = 0;
% end

%% Value function iteration

if verbose>=1
    fprintf('Start VFI.. Please be patient! \n')
end

V0=zeros(zdim,1);
V1=zeros(zdim,1);

apol = zeros(zdim,1);    % policy function for a-hat (index z=1,...,zdim)*/
%zpol = zeros(zdim,sdim); % policy for z+ depending on z for different s */


for i=1:zdim  % initialization (two options) */
    V0(i,1) = util(z_grid(i,1),gam)/(1-beta);
    %V0(i,1) = 0;
end

errt = 1;
iter = 1;

while errt>eps && iter<=500
    
    % Loop over z index (current state)
    start_a = 1;
    for iz = 1:zdim
        z_val = z_grid(iz,1);
        
        % Loop over all possible choices for a_hat'
        max_val = -inf;
        max_index = 1;
        for ia = start_a:zdim
            a_val  = a_grid(ia,1); % a_hat'
            % Compute continuation value
            Vnext = 0;
            for is_next = 1:sdim % future shocks s'
                znext = (1+rr)*a_val+ww*s_grid(is_next)-rr*phi;
                Vinterp = myinterp1(z_grid,V0,znext,1);
                Vnext = Vnext+pi_s(is_next)*Vinterp;
            end
            cons = z_val-a_val;
            if cons>0
                obj = util(cons,gam)+beta*Vnext;  % objective function: RHS of Bellman
                if obj>max_val
                    max_val = obj;
                    max_index = ia;
                else
                    if concavity==1; break; end
                end
            else
                break
            end % IF negative cons check
        end %end a_hat loop
        apol(iz) = a_grid(max_index,1);
        V1(iz)   = max_val;
        start_a  = max_index;
    end % end z loop
    
    errt = max(abs(V1-V0));
    % Update the current guess
    V0 = V1;
    
    if verbose>=2
        fprintf('iter = %d, err = %f \n',iter,errt)
    end
    
    iter = iter+1;
end

cpol = z_grid-apol;

if verbose>=1
    fprintf('Value function iteration converged after %d iter \n',iter)
end

if do_plots==1
    figure
    plot(z_grid,V1)
    title('Value function')
    
    gp = 20;
    figure
    plot(z_grid(1:gp),z_grid(1:gp),'--','linewidth',2)
    hold on
    plot(z_grid(1:gp),apol(1:gp),'linewidth',2)
    hold on
    plot(z_grid(1:gp),cpol(1:gp),'linewidth',2)
    legend('45','A(z)','C(z)')
    xlabel('Available resources z')
    ylabel('Consumption and assets')
    title('Policy function ahat')
end

%save temp

%% Simulation

if verbose>=1
    fprintf('Start simulation \n')
end

TT = 100000; % length of the simulation
s_sim = ones(TT,1);  % time series for exo shock s */
z_sim     = zeros(TT,1); % time series of z */
a_sim     = zeros(TT,1); % a-hat */
c_sim     = zeros(TT,1); % time series of c
asset_sim = zeros(TT,1); % assets=a-hat - Phic */

% Initial conditions
a_sim(1,1) = a_grid(1);

rng('default')
shock_vec = rand(TT,1);

for t = 2:TT
    shock = shock_vec(t);  % uncertainty revelation
    s_ind = ceil(shock*sdim);
    s_sim(t,1) = s_grid(s_ind);
    z_sim(t,1) = (1+rr)*a_sim(t-1,1)+ww*s_grid(s_ind)-rr*phi;
    a_sim(t,1) = myinterp1(z_grid,apol,z_sim(t,1),1);
    c_sim(t,1) = z_sim(t,1)-a_sim(t,1);
    asset_sim(t,1) = a_sim(t,1)-phi;
end

% Frequency distribution of assets,
% obtained from aatim without the first 100 iterations

asset_sim = asset_sim(101:end);
z_sim     = z_sim(101:end);
s_sim     = s_sim(101:end);
c_sim     = c_sim(101:end);

if do_plots==1
    figure
    subplot(3,1,1)
    plot(s_sim(1:300),'linewidth',2)
    title('Productivity shock, s')
    subplot(3,1,2)
    plot(asset_sim(1:300),'linewidth',2)
    title('Assets, a')
    subplot(3,1,3)
    plot(c_sim(1:300),'linewidth',2)
    title('Consumption, c')
end

if verbose>=1
    fprintf('Simulation is done! \n')
end

% Average assets
Ea = mean(asset_sim);

%% Update interest rate

%ED = sqrt(((Kd-Ea)/Kd)^2); % to be minimized
ED = Kd-Ea;

fprintf('New r = %f \n',rr)
fprintf('ED    = %f \n',ED)

disp('====================================================================')

% Pack outputs
sol.apol = apol;
sol.cpol = cpol;
sol.V1   = V1;

sim.s_sim = s_sim;
sim.asset_sim = asset_sim;
sim.c_sim = c_sim;
sim.z_sim = z_sim;

end

%-------------------------------------------------------------------------%
%                      INTERNAL PROCEDURES                                %
%-------------------------------------------------------------------------%

function F = util(c,gam)
% CRRA utility function. Assign very large negative if c is not feasible
F = -1e10;

if c>0
    F = c^(1-gam)/(1-gam);
end

end

function F = capital_demand(r,delta,alpha)
% Gives capital demand as a function of the interest rate

F = (alpha/(r+delta))^(1/(1-alpha));

end


function F = wage(K,alpha)
% Wage implied by firm's first order condition

F = (1-alpha)*K^alpha;

end

