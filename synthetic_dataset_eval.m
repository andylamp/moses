function synthetic_dataset_eval(n, T, r, alpha, nSim)
%%SYNTHETIC_DATASET_EVAL: This function is responsible for performing   
% comparisons for a number of synthetic datasets against three streaming  
% algorithms to compute the r-truncated SVD out of an incoming sequence 
% of vectors
%
% Author Andreas Grammenos (ag926@cl.cam.ac.uk)
%
% Last touched date 06/06/2018
% 
% License: GPLv3
%  

%% Initialise

% scope-in the global variables
global pflag;
global use_fast_moses_only
global use_offline_svds

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
[MosesT, MosesError, MosesFroT, MosesFT, MosesFError, MosesFFrotT, ...
  PowerT, PowerError, PowerFroT, GrouseT, GrouseError, GrouseFroT, ...
  OfflineT, OfflineError, OfflineFroT, ...
  Sigma] = online_svds_synthetic(n, r, T, alpha);

% Frobenius norm error normalised per block
PEs = nan(nSim, length(PowerT));
PEs(1,:) = PowerError;
GRs = nan(nSim, length(GrouseT));
GRs(1,:) = GrouseError;
% moses fast
MFEs = nan(nSim, length(MosesFT));
MFEs(1,:) = MosesFError;


%  MSE Errors
o_froT = nan(nSim, 1);
p_froT = nan(nSim, 1);
g_froT = nan(nSim, 1);
mf_froT = nan(nSim, 1);



% assign the first values
o_froT(1) = OfflineFroT;
p_froT(1) = PowerFroT;
g_froT(1) = GrouseFroT;
mf_froT(1) = MosesFFrotT;

if use_fast_moses_only == 0
  % moses simple
  MEs = nan(nSim, length(MosesT));
  MEs(1,:) = MosesError;
  
  m_froT = nan(nSim, 1);
  m_froT(1) = MosesFroT;
end

if use_offline_svds == 1
  OffEs = nan(nSim, length(OfflineT));
  OffEs(1, :) = OfflineError;
end

%profile off
%profile viewer
%pause

% loop for the remaining simulation
for i = 2:nSim
    fprintf("\n ** Simulation number %d **\n", i);
    [~, MosesError, MosesFroT, ~, MosesFError, MosesFFrotT, ...
     ~, PowerError, PowerFroT, ~, GrouseError, GrouseFroT, ...
     ~, OfflineError, OfflineFroT, ~] = online_svds_synthetic(n, r, T, alpha);
    
    % Frobenius norm error normalised with k^{.5}B
    PEs(i, :) = PowerError;
    GRs(i, :) = GrouseError;
    MFEs(i,:) = MosesFError;
    
    % Final Frobenius norm normalised with T error
    o_froT(i) = OfflineFroT;
    p_froT(i) = PowerFroT;
    g_froT(i) = GrouseFroT;
    mf_froT(i) = MosesFFrotT;
    
    if use_fast_moses_only == 0    
      MEs(i, :) = MosesError;
      m_froT(i) = MosesFroT;
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
AvGR = mean(GRs);

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
plot(GrouseT, AvGR, 'LineWidth', 2);

% plot offline svds only if we have to
if use_offline_svds == 1
  plot(OfflineT, AvOffE, 'LineWidth', 2);
end

hold off;
title('Rescaled error of Y_r relative to Y');
xlabel('samples'); ylabel('error');

% full legend cells
legendCells = {'MOSES', 'MOSES_s', 'PM', 'GROUSE'}; 

% remove moses simple if we are only running fast
if use_fast_moses_only == 1
  idc = ismember(legendCells, {'MOSES_s'});
  legendCells = legendCells(~idc);
end

% add the offline cell to the legends
if use_offline_svds == 1
  legendCells{end+1} = 'Offline';
end

% finally set the legends
legend(legendCells);

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
plot(GrouseT, AvGR, 'LineWidth', 2);

% plot offline svds only if we have to
if use_offline_svds == 1
  plot(OfflineT, AvOffE, 'LineWidth', 2);
end

hold off;
title('Rescaled error of Y_r relative to Y');
xlabel('samples'); ylabel('error');

% full legend cells
legendCells = {'MOSES', 'MOSES_s', 'PM', 'GROUSE'}; 

% remove moses simple if we are only running fast
if use_fast_moses_only == 1
  idc = ismember(legendCells, {'MOSES_s'});
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

%% Display only MOSES vs Power error over time

fig = figure;
hold on;
plot(MosesFT, AvMFE, 'LineWidth', 2);
plot(PowerT, AvPE, 'LineWidth', 2);
hold off;
legend('MOSES', 'PM');
title('MOSES vs PM Comparison');
t = sprintf("synthetic_froerror_moses_vs_pm_n_%s_T_%sk_r_%s_alpha_%s_nsim_%s", ...
  num2str(n), strrep(num2str(T/1000), ".", "_"), num2str(r), ...
  strrep(num2str(alpha), ".", "_"), ...
  strrep(num2str(nSim), ".", "_"));
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
plot(g_froT);
hold off;
title(['Error of final Y_r vs real Y over ' num2str(nSim) ' sims']);
xlabel('simulation number'); ylabel('error');

% set the legends accordingly
if use_fast_moses_only == 0
  legend('Offline', 'MOSES', 'MOSES_s', 'PM', 'GROUSE');
else
  legend('Offline', 'MOSES', 'PM', 'GROUSE');
end

subplot(2,1,2)
plot(mf_froT, 'LineWidth', 2);
hold on
plot(p_froT, 'LineWidth', 2);
legend('MOSES', 'PM');
title('Power Method vs MOSES'); % moses fast
xlabel('iterations'); ylabel('error');
hold off;

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
legend('MOSES', 'PM');
title('Power Method vs MOSES'); % moses fast
xlabel('iterations'); ylabel('error');
hold off;

% output figure to file if printing is enabled
t = sprintf("synthetic_fro_over_t_error_single_n_%s_T_%sk_r_%s_alpha_%s_nsim_%s", ...
  num2str(n), strrep(num2str(T/1000), ".", "_"), num2str(r), ...
  strrep(num2str(alpha), ".", "_"), ...
  strrep(num2str(nSim), ".", "_"));
print_fig(fig, t);

end

