function [MosesT, MosesError, MosesFroT, MosesFT, MosesFError, MosesFFroT,... 
  PowerT, PowerError, PowerFroT, ...
  FDT, FDError, FDFroT, ... 
  FDRT, FDRError, FDRFroT, ...
  GrouseT, GrouseError, GrouseFroT, ...
  OfflineT, OfflineError, OfflineFroT, ...
  Sigma] = online_svds_synthetic(n, r, T, alpha)
%ONLINE_SVDS_SYNTHETIC This function is responsible for performing a simulation
% comparing the three online SVD methods as well as the offline (if enabled).
%
% Author: Andreas Grammenos (ag926@cl.cam.ac.uk)
%
% Last touched date: 30/12/2018
% 
% License: GPLv3
%

global use_fast_moses_only
global use_offline_svds
global use_fdr

%% Data generation

[Y, ~, Sigma] = synthetic_data_gen(n, T, 1, alpha);

%% Offline r-leading SVD

[Soff, Doff, Voff, cflag] = svds(Y, 2*r);
% reduce the components to the r-truncation
Soff = Soff(:, 1:r); Doff = Doff(1:r, 1:r); Voff = Voff(:, 1:r);
% reconstruct Yr
Yroff = Soff * Doff * Voff';

% check for convergence
if cflag ~= 0
    fprintf("\n ** Offline SVD values did not converge (ret was: %d)\n", cflag);
else
    fprintf("\n ** Offline values did SVD converge (ret was: %d)\n", cflag);
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

[MosesFT, MosesFError, ~, ~, ~, Yr_mof, ~] = moses_fast(Y, r); % , 0, 2, 1

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

% only calculate that, if we have fdr enabled
if use_fdr == 1
  % enable error calculation for fd
  no_err = 0;
  % seed alpha
  a_seed = 0;
  % run the fdr
  [~, ~, FDRError, FDRT, Yr_fdr, ~] = fdr(Y', r, a_seed, no_err);
  % as with fd, since this is the transpose, revert it
  Yr_fdr = Yr_fdr';
else
  FDRT = NaN;
  FDRError = NaN;
  FDRFroT = NaN;
end

%% GROUSE method implemented from https://arxiv.org/pdf/1702.01005.pdf

[GrouseT, GrouseError, U_gr, V_gr, ~] = my_grouse(Y, r);

% reconstruct Yr based on grouse
Yr_gr = U_gr*V_gr';

%% Offline SVDS over all columns

if use_offline_svds == 1
  [OfflineT, OfflineError, ~] = incr_offline_svds(Y, r);
else
  OfflineT = 0;
  OfflineError = 0;
end

%% Frobenius norm normalised with T Error calculations

% Calculate the Frobenius normalised with T errors, namely: 
% n*Fro/T = n * Sum_{1}_{n} [ (Yr_i - Y_i)^2 ] / T

% min pad due to block alignment
min_pad = min([size(Yr_gr, 2), size(Yroff, 2), ...
  size(Yr_pm, 2), size(Yr_mof, 2), size(Yr_fd, 2)]);

% take in account moses simple, if needed
if use_fast_moses_only == 0
  min_pad = min([min_pad, size(Yr_mos, 2)]);
  MosesFroT = n*immse(Y(:, 1:min_pad), Yr_mos(:, 1:min_pad));
end

% take in account fdr, if needed
if use_fdr == 1
  min_pad = min([min_pad, size(Yr_fdr, 2)]);
  FDRFroT = n*immse(Y(:, 1:min_pad), Yr_fdr(:, 1:min_pad));
end

% calculate scaled mse (Fro err) & assign
OfflineFroT = n*immse(Y(:, 1:min_pad), Yroff(:, 1:min_pad));
PowerFroT = n*immse(Y(:, 1:min_pad), Yr_pm(:, 1:min_pad));
GrouseFroT = n*immse(Y(:, 1:min_pad), Yr_gr(:, 1:min_pad));
MosesFFroT = n*immse(Y(:, 1:min_pad), Yr_mof(:, 1:min_pad));
FDFroT = n*immse(Y(:, 1:min_pad), Yr_fd(:, 1:min_pad));

% Report them in a nice way
fprintf(" ** Final Frobenius norm over T Errors (Y vs YrHat)\n");
fprintf("\n\t -- Offline SVD: %d", OfflineFroT);
fprintf("\n\t -- Power Method: %d", PowerFroT);
fprintf("\n\t -- MOSES Faster: %d", MosesFFroT);
fprintf("\n\t -- FD: %d", FDFroT);
if use_fdr == 1
  fprintf("\n\t -- RFD: %d", FDRFroT);
end
if use_fast_moses_only == 0
  fprintf("\n\t -- MOSES: %d", MosesFroT);
end
fprintf("\n\t -- GROUSE: %d\n", GrouseFroT);
fprintf("\n");

end

