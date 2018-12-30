function [ Y, T, Sigma ] = synthetic_data_gen( n, T, lambda, alpha )
%SYNTHETIC_DATA_GEN function that generates T random vectors in R^{1 x n} 
% drawn from Power Law distribution with parameters `lambda` and `alpha`
%
% Author: Andreas Grammenos (ag926@cl.cam.ac.uk)
%
% Last touched date: 30/12/2018
% 
% License: GPLv3
% 

  % generate the singular spectrum
  dd = lambda*(1:n).^(-alpha);
  % generate Sigma
  Sigma = diag(dd);
  % random initialization of S basis
  S = orth(randn(n));
  % given S and Sigma generate the dataset (Y)
  Y = (S * Sigma * randn(n, T))/sqrt(T-1);

end

