function real_dataset_eval(path, r, desc)
%% This code implements a comparison for a real dataset against streaming  
% algorithms to compute rank-r truncated SVD of an incoming sequence of 
% vectors
%
% Author: Andreas Grammenos (ag926@cl.cam.ac.uk)
%
% Last touched date: 30/12/2018
% 
% License: GPLv3
%

%% Initialisation

% scope-in the global variables
global pflag;
global use_fast_moses_only
global use_fdr

% check if we have a description
if nargin < 3
  desc = "unknown";
  fprintf("\n\t ** Running real dataset evaluation\n");
else
  fprintf("\n\t ** Running real dataset evaluation using: %s dataset\n", desc);
end

% print iteration info
fprintf("\n\tTarget rank %d", r);
fprintf("\n\tPrint flag is: %d\n", pflag);

% load the dataset
Y = load(path);

% if number of columns in the default setting is larger than number of rows
% we use the transpose of Y
[numc, numr] = size(Y);
if numc > numr
  Y = Y';
end
% perform alignment
block_pad = 50;    % round up to the nearest block that is multiple of this
[rows, cols] = size(Y);
pad = mod(cols, block_pad);
Y = Y(:, 1: (cols-pad));

% update the column number
cols = size(Y, 2);

% center the dataset by using Y = Y - ((Y*(vec_ones*vec_ones'))./cols)
vec_ones = ones(cols, 1);
Y = Y - ((Y*(vec_ones*vec_ones'))./cols);


%% Run the comparison

[MosesT, MosesError, ~, ... 
  MosesFT, MosesFError, ~, ...
  PowerT, PowerError, ~,...
  FDT, FDError, ~, ...
  FDRT, FDRError, ~, ...
  GrouseT, GrouseError, ~, ~] = online_svds_real(Y, r);
  
%% Display error relative to T

fig = figure;
semilogy(MosesFT, MosesFError, 'LineWidth', 2);
hold on
if use_fast_moses_only == 0
  semilogy(MosesT, MosesError, 'LineWidth', 1);
end
semilogy(PowerT, PowerError, 'LineWidth', 2);
semilogy(FDT, FDError, 'LineWidth', 2);
if use_fdr == 1
  semilogy(FDRT, FDRError, 'LineWidth', 2);
end
semilogy(GrouseT, GrouseError, 'LineWidth', 2);
hold off;
caption = sprintf('Fro error over T of Y_r relative to Y for %s data', desc);
title(caption);
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

% finally set the legends
legend(legendCells, 'location', 'best');

% output figure to file if printing is enabled
t = sprintf("real_froerror_n_%s_r_%s_%s_dataset", ...
  num2str(rows), num2str(r), desc);
print_fig(fig, t);

end