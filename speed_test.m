function [] = speed_test(n, r, alpha, trials, T)
%SPEED_TEST performance testing of Moses, PM, and GROUSE
% performance evaluation based on the following criteria:
%
% n: the ambient dimension
% r: r-truncation target
% alpha: parameter for the power law distribution
% trials: the number of performance runs for each parameter tuple.
%
% We define as a parameter tuple the combination of:
%   ptuple = (n, [r, alpha]), where [r, alpha] is fixed and
%   n is changeable.
%
% Author Andreas Grammenos (ag926@cl.cam.ac.uk)
%
% Last touched date 06/06/2018
% 
% License: GPLv3
%  

%% Initialisation

% sanity check arguments
if size(n, 2) < 1 || 1 ~= size(r, 2) || 1 ~= size(alpha, 2) || trials <= 0
  fprintf("\n ** ERR: n must be an array, all others singular values, and trials must be > 0\n");
  return
end

% check for valid T (max time)
if nargin < 6
  % default max time is 10k  columns in R^(1 x n)
  T = 10000;
else
  if T < 2000
    fprintf("\n ** ERR: T must be at least 2k, reverting to default 10k\n");
    T = 10000;
  end
end

%% Run the trials

% raise the no error flag for speed tests
no_err = 1;
% set to the default floor multiplier
floor_mul = 2;
% number of parameter tuples
param_num = size(n, 2);

% preallocate time arrays
mos_ta = zeros(param_num, trials);
pm_ta = zeros(param_num, trials);
gr_ta = zeros(param_num, trials);

% now run the speed test
for i = 1:param_num
  fprintf("\n ** Starting running %d trials for parameter tuple (n=%d, r=%d ...)...\n", trials, n(i), r);
  % run the number of specified trials
  for j = 1:trials
    fprintf("\n -- Trial no: %d...\n", j);
    % first, generate the synthetic data for that trial
    Y = synthetic_data_gen(n(i), T, 1, alpha);
    % secondly, execute the timed runs
    [~, ~, ~, ~, ~, ~, mos_ta(i, j)] = moses_fast(Y, r, 2*r, floor_mul, no_err);
    [~, ~, ~, ~, pm_ta(i, j)] = mitliag_pm(Y, r, 2*n(i), floor_mul, no_err);
    [~, ~, ~, ~, gr_ta(i, j)] = my_grouse(Y, r, no_err);
  end
  fprintf("\nFinished running %d trials...\n", j);
end

%% Plot the results
mtr_mos = mean(mos_ta, 2);
mtr_pm = mean(pm_ta, 2);
mtr_gr = mean(gr_ta, 2);

fig = figure;
plot(mtr_mos, '-o', 'LineWidth', 2);
hold on
plot(mtr_pm, '-*', 'LineWidth', 2);
plot(mtr_gr, '-+', 'LineWidth', 2);
hold off
legend('MOSES', 'PM', 'GROUSE');
xticks(1:1:size(n, 2));
xticklabels(num2cell(n));
xlabel("ambient dimension (n)"); 
ylabel("average per trial time (sec)");
cap = sprintf("Time to compute Yr for r=%d", r);
title(cap);

% print figure, if needed
t = sprintf("speedtest_T_%sk_kr_%d_alpha_%d_trials_%d", ...
  strrep(num2str(T/1000), ".", "_"), r, alpha, trials);
print_fig(fig, t);

end

