function var = chop(var, tol)
%
% rounds the decimals var values up to the specified tol
%
% i.e. a = chop(1.23456789, 1e-3)
%      a = 1.235
%  
%
var = round(var./tol).*tol;
end