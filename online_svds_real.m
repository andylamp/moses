function [MosesT, MosesError, MosesFroT, MosesFT, MosesFError, MosesFFroT,... 
  PowerT, PowerError, PowerFroT, GrouseT, GrouseError, GrouseFroT, ...
  OfflineFroT] = online_svds_real(Y, r)
%%ONLINE_SVDS_REAL This function is responsible for performing a comparison
% of the three online SVD methods against the real-datasets.
%
% Author Andreas Grammenos (ag926@cl.cam.ac.uk)
%
% Last touched date 06/06/2018
% 
% License: GPLv3
% 

%% Initialisation

% scope-in the global variables
global use_fast_moses_only


%% Offline SVD of Y
n = size(Y, 1);

[Soff, Doff, Voff, cflag] = svds(Y, n); 
% reduce the components to the r-truncation
Soff = Soff(:, 1:r); Doff = Doff(1:r, 1:r); Voff = Voff(:, 1:r);
% get the offline svds reconstruction
Yroff = Soff * Doff * Voff';

if cflag ~= 0
    fprintf("\n\t ** Offline SVD values did not converge (ret was: %d)\n", cflag);
else
    fprintf("\n\t ** Offline values did SVD converge (ret was: %d)\n", cflag);
end

%% MOSES Simple (Our algorithm) from https://arxiv.org/pdf/1806.01304.pdf

if use_fast_moses_only == 0
  [MosesT, MosesError, ~, Yr_mos, ~] = moses_simple(Y, r);
else
  MosesT = NaN;
  MosesError = NaN;
  MosesFroT = NaN;
end

%% Moses Fast (Our algorithm, faster) from https://arxiv.org/pdf/1806.01304.pdf

[MosesFT, MosesFError, ~, ~, ~, Yr_mof, ~] = moses_fast(Y, r);

%% Power method implemented from https://arxiv.org/pdf/1307.0032.pdf

[PowerT, PowerError, ~, Yr_pm, ~] = mitliag_pm(Y, r);

%% GROUSE method implemented from https://arxiv.org/pdf/1702.01005.pdf

[GrouseT, GrouseError, U_gr, V_gr, ~] = my_grouse(Y, r);

% expand U_gr*V_gr' to get the Yr_gr
Yr_gr = U_gr*V_gr';

%% MSE error calculations

% Calculate the Frobenius normalised with T errors, namely: 
% n*Fro/T = n * Sum_{1}_{n} [ (Yr_i - Y_i)^2 ] / T

% find min offset
if use_fast_moses_only == 0
min_pad = min([size(Yroff, 2), size(Yr_pm, 2), ...
  size(Yr_mos, 2), size(Yr_mof, 2) size(Yr_gr, 2)]);
else
min_pad = min([size(Yroff, 2), size(Yr_pm, 2), ...
  size(Yr_mof, 2) size(Yr_gr, 2)]);
end

% extract the min pad version of Y
Y_aligned = Y(:, 1:min_pad);

% calculate
of_err = n*immse(Y_aligned, Yroff(:, 1:min_pad));
pm_err = n*immse(Y_aligned, Yr_pm(:, 1:min_pad));
gr_err = n*immse(Y_aligned, Yr_gr(:, 1:min_pad));
mof_err = n*immse(Y_aligned, Yr_mof(:, 1:min_pad));

% assign
OfflineFroT = of_err;
PowerFroT = pm_err;
GrouseFroT = gr_err;
MosesFFroT = mof_err;

if use_fast_moses_only == 0
  mo_err = n*immse(Y_aligned, Yr_mos(:, 1:min_pad));
  MosesFroT = mo_err;
end

% Report them in a nice way
fprintf(" ** Final Frobenius norm over T Errors (Y vs YrHat)\n");
fprintf("\n\t -- Offline SVD: %d", of_err);
fprintf("\n\t -- Power Method: %d", pm_err);
fprintf("\n\t -- MOSES Fast: %d", mof_err);
if use_fast_moses_only == 0
  fprintf("\n\t -- MOSES: %d", mo_err);
end
fprintf("\n\t -- GROUSE: %d\n", gr_err);
fprintf("\n");

end