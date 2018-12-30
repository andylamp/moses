function synthetic_dataset_eval(n, T, r, alpha, nSim)
%%SYNTHETIC_DATASET_EVAL: This function is responsible for performing   
% comparisons for a number of synthetic datasets against three streaming  
% algorithms to compute the r-truncated SVD out of an incoming sequence 
% of vectors
%
% Author: Andreas Grammenos (ag926@cl.cam.ac.uk)
%
% Last touched date: 30/12/2018
% 
% License: GPLv3
%

%% Initialise

% scope-in the global variables
global pflag;
global use_fast_moses_only
global use_offline_svds
global use_fdr

% sanity checks
if n > T
  fprintf("\n !! Error ambient dim (n) must be lower than T\n");
  return
end


%% Print iteration info

fprintf("\n ** Running evaluation with parameters:\n");
fprintf("\n\tPower law alpha=%d", alpha);
fprintf("\n\tAmbient dim: %d", n);
fprintf("\n\tTime (in number of Columns): %d", T);
fprintf("\n\tTarget rank: %d", r);
fprintf("\n\tPrint flag is: %d\n", pflag);

%% Run the simulation

%profile on

% run the initial simulation
fprintf("\n ** Simulation number 1 **\n");
 [MosesT, MosesError, MosesFroT, ...
   MosesFT, MosesFError, MosesFFrotT, ...
   PowerT, PowerError, PowerFroT, ...
   FDT, FDError, FDFroT, ...
   FDRT, FDRError, FDRFroT, ...
   GrouseT, GrouseError, GrouseFroT, ...
   OfflineT, OfflineError, OfflineFroT, ...
  Sigma] = online_svds_synthetic(n, r, T, alpha);

% Frobenius norm error normalised per block

% Power Method
PEs = nan(nSim, length(PowerT));
PEs(1, :) = PowerError;
% GROUSE
GREs = nan(nSim, length(GrouseT));
GREs(1, :) = GrouseError;
% Moses fast
MFEs = nan(nSim, length(MosesFT));
MFEs(1, :) = MosesFError;
% Frequent Directions
FDEs = nan(nSim, length(FDT));
FDEs(1, :) = FDError;
% Robust Frequent Directions
FDREs = nan(nSim, length(FDRT));
FDREs(1, :) = FDRError;

%  MSE Errors
o_froT = nan(nSim, 1);
p_froT = nan(nSim, 1);
g_froT = nan(nSim, 1);
mf_froT = nan(nSim, 1);
fd_froT = nan(nSim, 1);


% assign the first values
o_froT(1) = OfflineFroT;
p_froT(1) = PowerFroT;
g_froT(1) = GrouseFroT;
mf_froT(1) = MosesFFrotT;
fd_froT(1) = FDFroT;

% Only use that if we have Moses simple
if use_fast_moses_only == 0
  MEs = nan(nSim, length(MosesT));
  MEs(1,:) = MosesError;
  
  m_froT = nan(nSim, 1);
  m_froT(1) = MosesFroT;
end

% Only use that if we have fdr enabled
if use_fdr == 1
  fdr_froT = nan(nSim, 1);
  fdr_froT(1) = FDRFroT;
end

if use_offline_svds == 1
  % Offline SVD error
  OffEs = nan(nSim, length(OfflineT));
  OffEs(1, :) = OfflineError;
end

%profile off
%profile viewer
%pause

% loop for the remaining simulation
for i = 2:nSim
    fprintf("\n ** Simulation number %d **\n", i);
    [~, MosesError, MosesFroT, ...
     ~, MosesFError, MosesFFrotT, ...
     ~, PowerError, PowerFroT, ...
     ~, FDError, FDFroT, ...
     ~, FDRError, FDRFroT, ...
     ~, GrouseError, GrouseFroT, ...
     ~, OfflineError, OfflineFroT, ~] = online_svds_synthetic(n, r, T, alpha);
    
    % Frobenius norm error normalised with k^{.5}B
    PEs(i, :) = PowerError;
    GREs(i, :) = GrouseError;
    MFEs(i, :) = MosesFError;
    FDEs(i, :) = FDError;
    
    % Final Frobenius norm normalised with T error
    o_froT(i) = OfflineFroT;
    p_froT(i) = PowerFroT;
    g_froT(i) = GrouseFroT;
    mf_froT(i) = MosesFFrotT;
    fd_froT(i) = FDFroT;
    
    if use_fast_moses_only == 0    
      MEs(i, :) = MosesError;
      m_froT(i) = MosesFroT;
    end
    
    if use_fdr == 1
      FDREs(i, :) = FDRError;
      fdr_froT(i) = FDRFroT;
    end
    
    if use_offline_svds == 1
      OffEs(i, :) = OfflineError;
    end
end

% calculate relative metrics

% moses fast averaged metrics
AvMFE = mean(MFEs);

% power method averaged metrics
AvPE = mean(PEs);

% grouse metrics
AvGR = mean(GREs);

% fd metrics
AvFD = mean(FDEs);

% fdr metrics
AvFDR = mean(FDREs);

% moses simple averaged metrics  
if use_fast_moses_only == 0
  AvME = mean(MEs);
end

% offline averaged metrics
if use_offline_svds == 1
  AvOffE = mean(OffEs);
end

%% Display scree plot of population distribution

fig = figure;
subplot(2,1,1)
plot(1:n, sqrt(diag(Sigma).^2));
title(['Population scree plot (alpha: ' num2str(alpha) ')'])
xlabel('rank'); ylabel('singular value')



%% Display error relative to T

subplot(2,1,2)
plot(MosesFT, AvMFE, 'LineWidth', 2);
hold on

% plot moses simple only if we have to
if use_fast_moses_only == 0
  plot(MosesT, AvME, 'LineWidth', 1);
end
plot(PowerT, AvPE, 'LineWidth', 2);
plot(FDT, AvFD, 'LineWidth', 2);

% plot fdr only if we have to
if use_fdr == 1
  plot(FDRT, AvFDR, 'LineWidth', 2);
end

plot(GrouseT, AvGR, 'LineWidth', 2);

% plot offline svds only if we have to
if use_offline_svds == 1
  plot(OfflineT, AvOffE, 'LineWidth', 2);
end

hold off;
title('Rescaled error of Y_r relative to Y');
xlabel('samples'); ylabel('error');

% full legend cells
legendCells = {'MOSES', 'MOSES_s', 'PM', 'FD', 'FDR', 'GROUSE'}; 

% remove moses simple if we are only running fast
if use_fast_moses_only == 1
  idc = ismember(legendCells, {'MOSES_s'});
  legendCells = legendCells(~idc);
end

% remove fdr if need be
if use_fdr == 0
  idc = ismember(legendCells, {'FDR'});
  legendCells = legendCells(~idc);
end

% add the offline cell to the legends
if use_offline_svds == 1
  legendCells{end+1} = 'Offline';
end

% finally set the legends
legend(legendCells, 'Location', 'best');

% set the figure limits
xlim([1, T])
% set the axis labels correctly
xlabel('samples'); ylabel('error');

% output figure to file if printing is enabled
t = sprintf("synthetic_froerror_n_%s_r_%s_alpha_%s_nsim_%s", ...
  num2str(n), num2str(r), ...
  strrep(num2str(alpha), ".", "_"), ...
  strrep(num2str(nSim), ".", "_"));
print_fig(fig, t);

%% Display the error relative to T as a singular plot

fig = figure;
plot(MosesFT, AvMFE, 'LineWidth', 2);
hold on

% plot moses simple only if we have to
if use_fast_moses_only == 0
  plot(MosesT, AvME, 'LineWidth', 1);
end
plot(PowerT, AvPE, 'LineWidth', 2);
plot(FDT, AvFD, 'LineWidth', 2);

if use_fdr == 1
  plot(FDRT, AvFDR, 'LineWidth', 2);
end
plot(GrouseT, AvGR, 'LineWidth', 2);

% plot offline svds only if we have to
if use_offline_svds == 1
  plot(OfflineT, AvOffE, 'LineWidth', 2);
end

hold off;
title('Rescaled error of Y_r relative to Y');
xlabel('samples'); ylabel('error');

% full legend cells
legendCells = {'MOSES', 'MOSES_s', 'PM', 'FD', 'FDR', 'GROUSE'}; 

% remove moses simple if we are only running fast
if use_fast_moses_only == 1
  idc = ismember(legendCells, {'MOSES_s'});
  legendCells = legendCells(~idc);
end

% remove fdr if need be
if use_fdr == 0
  idc = ismember(legendCells, {'FDR'});
  legendCells = legendCells(~idc);
end

% add the offline cell to the legends
if use_offline_svds == 1
  legendCells{end+1} = 'Offline';
end

% plot the correct legends
legend(legendCells, 'Location', 'SouthEast');

% set the figure limits
xlim([1, T])
% set the axis labels correctly
xlabel('samples'); ylabel('error');
% output figure to file if printing is enabled
t = sprintf("synthetic_froerror_noscree_n_%s_r_%s_alpha_%s_nsim_%s", ...
  num2str(n), num2str(r), ...
  strrep(num2str(alpha), ".", "_"), ...
  strrep(num2str(nSim), ".", "_"));
print_fig(fig, t);

%% Display only MOSES vs Power vs FD error over time

fig = figure;
hold on;
plot(MosesFT, AvMFE, 'LineWidth', 2);
plot(PowerT, AvPE, 'LineWidth', 2);
plot(FDT, AvFD, 'LineWidth', 2);

% plot fdr if enabled
if use_fdr == 1
  plot(FDRT, AvFDR, 'LineWidth', 2);
end
hold off;

if use_fdr == 1
  legend('MOSES', 'PM', 'FD', 'FDR');
  title('MOSES vs PM vs FD vs FDR Comparison');
  t = sprintf("synthetic_froerror_moses_vs_pm_vs_fd_vs_fdr_n_%s_T_%sk_r_%s_alpha_%s_nsim_%s", ...
    num2str(n), strrep(num2str(T/1000), ".", "_"), num2str(r), ...
    strrep(num2str(alpha), ".", "_"), ...
    strrep(num2str(nSim), ".", "_"));
else
  legend('MOSES', 'PM', 'FD');
  title('MOSES vs PM vs FD Comparison');
  t = sprintf("synthetic_froerror_moses_vs_pm_vs_fd_n_%s_T_%sk_r_%s_alpha_%s_nsim_%s", ...
    num2str(n), strrep(num2str(T/1000), ".", "_"), num2str(r), ...
    strrep(num2str(alpha), ".", "_"), ...
    strrep(num2str(nSim), ".", "_"));
end
print_fig(fig, t);


%% Display error relative to the Frobenius norm over T of final Y_r against Y

fig = figure;
subplot(2,1,1)
hold on
plot(o_froT);
plot(mf_froT);
if use_fast_moses_only == 0
  plot(m_froT);
end
plot(p_froT);
plot(fd_froT);

% plot fdr only we have to
if use_fdr == 1
  plot(fdr_froT);
end
plot(g_froT);
hold off;
title(['Error of final Y_r vs real Y over ' num2str(nSim) ' sims']);
xlabel('simulation number'); ylabel('error');

% full legend cells
legendCells = {'MOSES', 'MOSES_s', 'PM', 'FD', 'FDR', 'GROUSE'}; 

% remove moses simple if we are only running fast
if use_fast_moses_only == 1
  idc = ismember(legendCells, {'MOSES_s'});
  legendCells = legendCells(~idc);
end

% remove fdr if need be
if use_fdr == 0
  idc = ismember(legendCells, {'FDR'});
  legendCells = legendCells(~idc);
end

% finally set the legend cells
legend(legendCells);

subplot(2,1,2)
plot(mf_froT, 'LineWidth', 2);
hold on
plot(p_froT, 'LineWidth', 2);
plot(fd_froT, 'LineWidth', 2);

% plot fdr only we have to
if use_fdr == 1
  plot(fdr_froT, 'LineWidth', 2);
end

hold off;

% full legend cells
legendCells = {'MOSES', 'PM', 'FD', 'FDR'}; 

% remove fdr if need be
if use_fdr == 0
  idc = ismember(legendCells, {'FDR'});
  legendCells = legendCells(~idc);
  title('Power Method vs MOSES vs FD');
else
  title('Power Method vs MOSES vs FD vs FDR');
end

legend(legendCells);
xlabel('iterations'); ylabel('error');

% output figure to file if printing is enabled
t = sprintf("synthetic_fro_over_t_error_n_%s_T_%sk_r_%s_alpha_%s_nsim_%s", ...
  num2str(n), strrep(num2str(T/1000), ".", "_"), num2str(r), ...
  strrep(num2str(alpha), ".", "_"), ...
  strrep(num2str(nSim), ".", "_"));
print_fig(fig, t);

%% Display error relative to the Frobenius norm over T of final Y_r against Y 
% for only MOSES and PM

fig = figure;
plot(mf_froT, 'LineWidth', 2);
hold on
plot(p_froT, 'LineWidth', 2);
plot(fd_froT, 'LineWidth', 2);

% plot fdr only we have to
if use_fdr == 1
  plot(fdr_froT, 'LineWidth', 2);
end
hold off;

% full legend cells
legendCells = {'MOSES', 'PM', 'FD', 'FDR'}; 

% remove fdr if need be
if use_fdr == 0
  idc = ismember(legendCells, {'FDR'});
  legendCells = legendCells(~idc);
  title('MOSES vs Power Method vs FD');
else
  title('MOSES vs Power Method vs FD vs FDR');
end

legend(legendCells);
xlabel('iterations'); ylabel('error');

% output figure to file if printing is enabled
t = sprintf("synthetic_fro_over_t_error_single_n_%s_T_%sk_r_%s_alpha_%s_nsim_%s", ...
  num2str(n), strrep(num2str(T/1000), ".", "_"), num2str(r), ...
  strrep(num2str(alpha), ".", "_"), ...
  strrep(num2str(nSim), ".", "_"));
print_fig(fig, t);

end

