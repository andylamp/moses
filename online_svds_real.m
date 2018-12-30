function [MosesT, MosesError, MosesFroT, ...
  MosesFT, MosesFError, MosesFFroT,... 
  PowerT, PowerError, PowerFroT, ...
  FDT, FDError, FDFroT, ... 
  FDRT, FDRError, FDRFroT, ... 
  GrouseT, GrouseError, GrouseFroT, ...
  OfflineFroT] = online_svds_real(Y, r)
%%ONLINE_SVDS_REAL This function is responsible for performing a comparison
% of the three online SVD methods against the real-datasets.
%
% Author: Andreas Grammenos (ag926@cl.cam.ac.uk)
%
% Last touched date: 30/12/2018
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

%% Frequent Directions method as seen in https://arxiv.org/abs/1501.01711.pdf

% enable error calculation for fd
no_err = 0;
% run the fd
[~, FDError, FDT, Yr_fd, ~] = fd(Y', r, no_err);
% since this is the transpose, revert it
Yr_fd = Yr_fd';

%% Robust Frequent Direction method as seen in https://arxiv.org/pdf/1705.05067

% enable error calculation for fd
no_err = 0;
% seed alpha
a_seed = 0;
% run the fdr
[~, ~, FDRError, FDRT, Yr_fdr, ~] = fdr(Y', r, a_seed, no_err);
% as with fd, since this is the transpose, revert it
Yr_fdr = Yr_fdr';

%% GROUSE method implemented from https://arxiv.org/pdf/1702.01005.pdf

[GrouseT, GrouseError, U_gr, V_gr, ~] = my_grouse(Y, r);

% expand U_gr*V_gr' to get the Yr_gr
Yr_gr = U_gr*V_gr';

%% MSE error calculations

% Calculate the Frobenius normalised with T errors, namely: 
% n*Fro/T = n * Sum_{1}_{n} [ (Yr_i - Y_i)^2 ] / T

% find min offset
min_pad = min([size(Yroff, 2), size(Yr_pm, 2), ...
  size(Yr_mof, 2) size(Yr_gr, 2), size(Yr_fd, 2), ...
  size(Yr_fdr, 2)]);

% take in account moses simple, if needed
if use_fast_moses_only == 0
    min_pad = min([min_pad, size(Yr_mos, 2)]);
    MosesFroT = n*immse(Y(:, 1:min_pad), Yr_mos(:, 1:min_pad));
end

% extract the min pad version of Y
Y_aligned = Y(:, 1:min_pad);

% calculate scaled mse (Fro err) & assign
OfflineFroT = n*immse(Y_aligned, Yroff(:, 1:min_pad));
PowerFroT = n*immse(Y_aligned, Yr_pm(:, 1:min_pad));
GrouseFroT = n*immse(Y_aligned, Yr_gr(:, 1:min_pad));
MosesFFroT = n*immse(Y_aligned, Yr_mof(:, 1:min_pad));
FDFroT = n*immse(Y_aligned, Yr_fd(:, 1:min_pad));
FDRFroT = n*immse(Y_aligned, Yr_fdr(:, 1:min_pad));

% Report them in a nice way
fprintf(" ** Final Frobenius norm over T Errors (Y vs YrHat)\n");
fprintf("\n\t -- Offline SVD: %d", OfflineFroT);
fprintf("\n\t -- Power Method: %d", PowerFroT);
fprintf("\n\t -- MOSES Fast: %d", MosesFFroT);
fprintf("\n\t -- FD: %d", FDFroT);
fprintf("\n\t -- FDR: %d", FDRFroT);
if use_fast_moses_only == 0
  fprintf("\n\t -- MOSES: %d", MosesFroT);
end
fprintf("\n\t -- GROUSE: %d\n", GrouseFroT);
fprintf("\n");

end