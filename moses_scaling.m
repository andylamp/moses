function [] = moses_scaling(n_arr, r_arr, bmul_arr, T, alpha, trials)
%MOSES_SCALING MOSES_fast scaling experiment for different n, r, blk
%permutations.
%
% Author: Andreas Grammenos (ag926@cl.cam.ac.uk)
%
% Last touched date: 30/12/2018
% 
% License: GPLv3
% 

%% Initialise

% scope in global variables
global run_full_scaling
global run_exp1
global run_exp2
global run_exp3

% check if we have T
if nargin < 4
  T = 3000;  % default T
end

% check if we have alpha
if nargin < 5
  alpha = 1; % default alpha
end

% check if we have a trials value
if nargin < 6
  trials = 10;
end

%% Generic run for all tuples

% Run for all parameter tuples, in this instance we define a parameter
% tuple as p = (n, r, b) and each of the graphs are generated by fixing two
% parameters for any given plot while the third one changes.

if run_full_scaling == 1
  fprintf("\n ** Running generic scaling experiment **\n\n");

  % print the configuration parameters in a nice way
  fprintf("\n -- Provided Configuration Parameters\n");
  fprintf("\n\t T: %d, alpha %d", T, alpha);
  fprintf("\n\t Ambient dimensions: %s", num2str(n_arr));
  fprintf("\n\t Target ranks: %s", num2str(r_arr));
  fprintf("\n\t Block size multipliers: %s\n", num2str(bmul_arr));
  
  % preallocate fro_err
  fro_err = zeros(size(n_arr, 2), size(r_arr, 2), size(bmul_arr, 2));
  % run for all permutations -- EXPENSIVE!
  for i = 1:size(n_arr, 2)
    [Y, ~, ~] = synthetic_data_gen(n_arr(i), T, 1, alpha);
    for j = 1:size(r_arr, 2)
      fprintf("\n !! Running for n=%d, T=%dk, r=%d\n", ...
        n_arr(i), T/1000, r_arr(j));
      fig = figure;
      hold on
      for k = 1:size(bmul_arr, 2) 
        blk = bmul_arr(k)*r_arr(j);
        [mft, m_err, ~, ~, ~, Yr_mof, ~] = moses_fast(Y, r_arr(j), blk);
        % check if we need to pad due to block misalignment
        min_pad = size(Yr_mof, 2);
        % calculate the final error
        fro_err(i, j, k) = n_arr(i)*immse(Y(:, 1:min_pad), Yr_mof);
        % plot the error
        plot(mft, m_err, 'LineWidth', 2);
      end
      % construct the title
      cap = sprintf("moses_f error n=%d, r=%d", n_arr(i), r_arr(j));
      % assign the title
      title(cap);
      % construct the legends
      leg = cellstr(num2str(bmul_arr', 'b=%-dr'));
      % assign the legends
      legend(leg);    
      hold off
      % now print
      t = sprintf("moses_scaling_fro_over_t_error_n_%s_T_%sk_r_%s", ...
        num2str(n_arr(i)), strrep(num2str(T/1000), ".", "_"), ...
        num2str(r_arr(j)));
      print_fig(fig, t);
    end
    fprintf("\n ** Finished running generic scaling experiment **\n");
  end
else
  fprintf("\n ** Skipping running generic scaling experiment **\n");
end

%% Experiment one (r = 15, b = 2r for all n_arr)

if run_exp1 == 1
  
  fprintf("\n ** Running experiment 1 (fixed r, b -- variable n)\n");

  r = 15;                             % r-recovery
  m_err_a = zeros(1, size(n_arr, 2)); % holds the final errors
  fig = figure;
  hold on;
  for i = 1:size(n_arr, 2)
    fprintf("\n !! Running for n=%d\n", n_arr(i));
    for j = 1:trials
      fprintf("\n !! Trial %d out of %d\n", j, trials);
      [Y, ~, ~] = synthetic_data_gen(n_arr(i), T, 1, alpha);
      [mft, m_err, ~, ~, ~, ~, ~] = moses_fast(Y, r);
      if j == 1
        % plot the error
        plot(mft, m_err, 'LineWidth', 2);
      end
      % find the last non-NaN location of the error array
      l_idc = find(sum(~isnan(m_err),1) > 0, 1 , 'last');
      % append the final error we have on record
      m_err_a(i) = m_err(l_idc);
    end
  end
  hold off;
  % construct the title
  cap = sprintf("MOSES fixed b=2r, r=%d variable n", r);
  % assign the title
  title(cap);
  % construct the legends
  leg = cellstr(num2str(n_arr', 'n=%-d'));
  % assign the legends
  legend(leg);    
  % assign axis labels
  ylabel("error"); xlabel("samples");   
  hold off

  % now print
  t = sprintf("moses_scaling_exp1_fro_over_t_error_fixed_n_r_var_n_T_%sk", ...
    strrep(num2str(T/1000), ".", "_"));
  print_fig(fig, t);

  % now plot the final error by itself
  fig = figure;
  plot(m_err_a./trials, '-*', 'LineWidth', 2);
  % construct the title
  cap = sprintf("MOSES fixed b=2r, r=%d variable n final error", r);
  title(cap);
  % put the correct axis labels
  xlabel("ambient dimension (n)"); ylabel("error");
  % set the correct x-axis ticks as well as their labels
  xticks(1:1:size(n_arr, 2));
  xticklabels(num2cell(n_arr));
  % increase the font size of the figure
  set(gca,'fontsize', 18);

  % now print
  t = sprintf("moses_scaling_exp1_final_error_fixed_n_r_var_n_T_%sk", ...
    strrep(num2str(T/1000), ".", "_"));
  print_fig(fig, t);

fprintf("\n ** Finished running experiment 1\n");

else
  fprintf("\n ** Skipping running experiment 1\n");
end

%% Experiment two (r = 15, n = max(n_arr), variable b)

if run_exp2 == 1

  fprintf("\n ** Running experiment 2 (fixed r, n -- variable b)\n");

  r = 15;                                 % r-recovery
  n = max(n_arr);                         % ambient dim
  m_err_a = zeros(1, size(bmul_arr, 2));  % holds the final errors
  fig = figure;
  hold on;
  for i = 1:size(bmul_arr, 2)
    fprintf("\n !! Running for b=%dr\n", bmul_arr(i));
    for j = 1:trials
      fprintf("\n\t ** Trial %d out of %d\n", j, trials);
      [Y, ~, ~] = synthetic_data_gen(n, T, 1, alpha);
      [mft, m_err, ~, ~, ~, ~, ~] = moses_fast(Y, r, bmul_arr(i)*r);
      % plot the error
      if j == 1
        plot(mft, m_err, 'LineWidth', 2);
      end
      % find the last non-NaN location of the error array
      l_idc = find(sum(~isnan(m_err),1) > 0, 1 , 'last');
      % append the final error we have on record
      m_err_a(i) = m_err_a(i) + m_err(l_idc);
    end
  end
  hold off;
  % construct the title
  cap = sprintf("MOSES fixed n=%d, r=%d variable b", n, r);
  % assign the title
  title(cap);
  % construct the legends
  leg = cellstr(num2str(bmul_arr', 'b=%-dr'));
  % assign the legends
  legend(leg);    
  % set the legends
  xlabel('samples'); ylabel('error');
  hold off
  % now print
  t = sprintf("moses_scaling_exp2_fro_over_t_error_fixed_n_r_var_b_T_%sk", ...
    strrep(num2str(T/1000), ".", "_"));
  print_fig(fig, t);

  % now plot the final error by itself
  fig = figure;
  plot(m_err_a./trials, '-*', 'LineWidth', 2);
  % construct the title
  cap = sprintf("MOSES fixed n=%d, r=%d variable b final error", n, r);
  title(cap);
  % put the correct axis labels
  xlabel("block size (b)"); ylabel("error");
  % set the correct x-axis ticks as well as their labels
  xticks(1:2:size(bmul_arr, 2));
  ticks = num2str(bmul_arr(3:2:end)', '%-dr');
  ftick = blanks(size(ticks, 2));
  ftick(end) = 'r';
  ticks = [ftick; ticks];
  xticklabels(cellstr(ticks));
  xlim([1, size(bmul_arr, 2)]);
  % increase the font size of the figure
  set(gca,'fontsize', 18);

  % now print
  t = sprintf("moses_scaling_exp2_final_error_fixed_n_r_var_b_T_%sk", ...
    strrep(num2str(T/1000), ".", "_"));
  print_fig(fig, t);

  fprintf("\n ** Finished running experiment 2\n");
else
  
  fprintf("\n ** Skipping running experiment 2\n");
end

%% Experiment three (n = max(n_arr), b = 2r, variable r)

if run_exp3 == 1

  fprintf("\n ** Running experiment 3 (fixed n, b -- variable r)\n");

  n = max(n_arr);                       % ambient dim
  m_err_a = zeros(1, size(r_arr, 2));   % holds the final errors
  fig = figure;
  hold on;
  for i = 1:size(r_arr, 2)
    fprintf("\n !! Running for r=%d\n", r_arr(i));
    for j = 1:trials
      fprintf("\n !! Trial %d out of %d\n", j, trials);
      [Y, ~, ~] = synthetic_data_gen(n, T, 1, alpha);
      [mft, m_err, ~, ~, ~, ~, ~] = moses_fast(Y, r_arr(i));
      if j == 1
        % plot the error
        plot(mft, m_err, 'LineWidth', 2);
      end
      % find the last non-NaN location of the error array
      l_idc = find(sum(~isnan(m_err),1) > 0, 1 , 'last');
      % append the final error we have on record
      m_err_a(i) = m_err(l_idc);
    end
  end
  hold off;
  % construct the title
  cap = sprintf("MOSES fixed n=%d, b=2r variable r", n);
  % assign the title
  title(cap);
  % construct the legends
  leg = cellstr(num2str(r_arr', 'r=%-d'));
  % assign the legends
  legend(leg); 
  % set the legends
  xlabel('samples'); ylabel('error');   
  hold off
  % now print
  t = sprintf("moses_scaling_exp3_fro_over_t_error_fixed_n_b_var_r_T_%sk", ...
    strrep(num2str(T/1000), ".", "_"));
  print_fig(fig, t);

  % now plot the final error by itself
  fig = figure;
  plot(m_err_a./trials, '-*', 'LineWidth', 2);
  % construct the title
  cap = sprintf("MOSES fixed n=%d, b=2r variable r final error", n);
  title(cap);
  % put the correct axis labels
  xlabel("rank (r)"); ylabel("error");
  % set the correct x-axis ticks as well as their labels
  xticks(1:1:size(r_arr, 2));
  xticklabels(num2cell(r_arr));
  % increase the font size of the figure
  set(gca,'fontsize', 18);

  % now print
  t = sprintf("moses_scaling_exp3_final_error_fixed_n_b_var_r_T_%sk", ...
    strrep(num2str(T/1000), ".", "_"));
  print_fig(fig, t);

  fprintf("\n ** Finished running experiment 3\n");
else
  fprintf("\n ** Skipped running experiment 3\n");
end

