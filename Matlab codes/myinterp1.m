function yi = myinterp1(x,y,xi,extrap)

% INPUTS:
% x and y are N*1 vectors
% x MUST be monotonically increasing
% xi is the query point and it MUST be a scalar
% extrap is a 0-1 integer
% To do extrapolation, set extrap=1.
% Always recommended to set extrap=1! 

% OUTPUT
% yi is the value of the interpolated function at query point xi 

% -------------------------------------------------------------------------
% Copyright © 2018 by Alessandro Di Nola. All rights 
% reserved.
% -------------------------------------------------------------------------

n = size(x,1);
if size(y,1)~=n
    error('myinterp1: x and y must be vectors with same length!')
end

% Find x(j)<= xi <x(j+1), for j=1,..,n-1

%% Using Paul Klein's locate
j = locate(x,xi);
% xi is between x(j) and xx(j+1)
%
% j = 0 and j = n means x is out of range!!!!!!!!!!!!!!

j = max(min(j,n-1),1);

slope = (y(j+1)-y(j))/(x(j+1)-x(j));
yi = y(j)+(xi-x(j))*slope;

% Give NaN to the values of 'yi' corresponding to 'xi' out of the range of 'x'
if extrap==0
    yi(xi<x(1) | xi>x(end)) = NaN;
end

end