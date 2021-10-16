function [Tran,s,p,arho,asigma]=markovapprox(rho,sigma,mu,m,N,disp_on_screen)

%Syntax: [Tran,s,p,arho,asigma]=markovappr(rho,sigma,m,N)
%
%This function pproximates a first-order autoregressive process
%with persistence rho and innovation standard deviation sigma with
%an N state Markov chain; m determines the width of discretized state
%space, Tauchen uses m=3, with ymax=m*vary,ymin=-m*vary, where ymax
%and ymin are the two boundary points, Tran is the transition matrix
%of the Markov chain, s is the discretized state space, p is the
%chain stationary distribution, arho is the theoretical first order
%autoregression coefficient for the Markov chain, asigma is the
%theoretical standard deviation for the Markov chain.
%
%Eva Carceles-Poveda 2003

if nargin<6
    disp_on_screen = 0;
end

%Discretize the state space
stvy = sqrt(sigma^2/(1-rho^2)); 	% standard deviation of y(t)
ystar = mu/(1.0d0-rho);             % expected value of y
ymax = m*stvy;                    	% upper boundary of state space
ymin = -ymax;                     	% lower boundary of state space
w = (ymax-ymin)/(N-1);    	     	% distance between points
s = ystar + (ymin:w:ymax);          % the discretized state space


% Calculate the transition matrix
Tran = zeros(N,N);
for j=1:N
    for k=2:N-1
        Tran(j,k)= normcdf(s(k)-rho*s(j)+w/2,mu,sigma)...
            - normcdf(s(k)-rho*s(j)-w/2,mu,sigma);
    end
    Tran(j,1) = normcdf(s(1)-rho*s(j)+w/2,mu,sigma);
    Tran(j,N) = 1 - normcdf(s(N)-rho*s(j)-w/2,mu,sigma);
end

% Check that Tran is well specified
sumRows = sum(Tran,2);
if max(abs( sumRows - 1 )) > 1e-6
    error('Probabilities do not sum to one');
end

% if sum(Tran') ~= ones(1,N)
%    str = find(Tran'-ones(1,N));  % find rows not adding up to one
%    disp('error in transition matrix');
%    disp(['rows ',num2str(str),' does not sum to one']);
% end


% Calculate the invariant distribution of Markov chain
Trans= Tran';
p = (1/N)*ones(N,1); % initial distribution of states
test = 1;
while test > 10^(-8)
    p1 = Trans*p;
    test=max(abs(p1-p));
    p = p1;
end

meanm = s*p;             		  % mean of invariant distribution of chain
varm = ((s-meanm).^2)*p;  		  % variance of invariant distribution of chain
midaut1 = (s-meanm)'*(s-meanm);   % cross product of deviation from mean of yt and yt-1
probmat = p*ones(1,N);     		  % each column is invariant distribution
midaut2 = Tran.*probmat.*midaut1; % product of the first two terms is joint distribution of (Yt-1,Yt)
autcov1 = sum(sum(midaut2));      %  first-order auto-covariance


if (disp_on_screen==1)
    
    % Display moments of chain
    disp('rho of original process v.s that of Markov chain')
    disp('')
    arho = autcov1/varm;            % theoretical rho
    disp([rho arho])
    
    disp('standard deviation of true process v.s that of Markov chain')
    disp('')
    asigma = sqrt(varm);
    disp([stvy asigma])
    
end

s = s';

end %end function <markovapprox>
