function [MosesT, MosesError, MosesFroT, MosesFT, MosesFError, MosesFFroT,... 
  PowerT, PowerError, PowerFroT, GrouseT, GrouseError, GrouseFroT, ...
  OfflineT, OfflineError, OfflineFroT, ...
  Sigma] = online_svds_synthetic(n, r, T, alpha)
%ONLINE_SVDS_SYNTHETIC This function is responsible for performing a simulation
% comparing the three online SVD methods as well as the offline (if enabled).
%
% Author Andreas Grammenos (ag926@cl.cam.ac.uk)
%
% Last touched date 06/06/2018
% 
% License: GPLv3
% 

global use_fast_moses_only
global use_offline_svds

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
   

%% MOSES Simple (Our algorithm) from https://arxiv.org/abs/1806.01304

if use_fast_moses_only == 0
  [MosesT, MosesError, ~, Yr_mos, ~] = moses_simple(Y, r);
else
  MosesT = NaN;
  MosesError = NaN;
  MosesFroT = NaN;
end

%% Moses Fast (Our algorithm, faster) from https://arxiv.org/abs/1806.01304

[MosesFT, MosesFError, ~, ~, ~, Yr_mof, ~] = moses_fast(Y, r);

%% Power method implemented from https://arxiv.org/pdf/1307.0032.pdf

[PowerT, PowerError, ~, Yr_pm, ~] = mitliag_pm(Y, r);

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
pad = min([size(Yr_gr, 2), size(Yroff, 2), ...
  size(Yr_pm, 2), size(Yr_mof, 2)]);

if use_fast_moses_only == 0
  % take in account moses simple, if needed
  pad = min([pad, size(Yr_mos, 2)]);
  mos_err = n*immse(Y(:, 1:pad), Yr_mos(:, 1:pad));
  MosesFroT = mos_err;
else
  
end

% calculate
of_err = n*immse(Y(:, 1:pad), Yroff(:, 1:pad));
pm_err = n*immse(Y(:, 1:pad), Yr_pm(:, 1:pad));
gr_err = n*immse(Y(:, 1:pad), Yr_gr(:, 1:pad));
mof_err = n*immse(Y(:, 1:pad), Yr_mof(:, 1:pad));

% assign
OfflineFroT = of_err;
PowerFroT = pm_err;
GrouseFroT = gr_err;
MosesFFroT = mof_err;

% Report them in a nice way
fprintf(" ** Final Frobenius norm over T Errors (Y vs YrHat)\n");
fprintf("\n\t -- Offline SVD: %d", of_err);
fprintf("\n\t -- Power Method: %d", pm_err);
fprintf("\n\t -- MOSES Faster: %d", mof_err);
if use_fast_moses_only == 0
  fprintf("\n\t -- MOSES: %d", mos_err);
end
fprintf("\n\t -- GROUSE: %d\n", gr_err);
fprintf("\n");

end

