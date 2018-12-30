%% Comparison script for Streaming r-truncated SVD (Moses, PM, FD, & GROUSE)
%
% Description: 
%   This code is supplied as additional material alongside our paper:
%   "MOSES: A Streaming Algorithm for Linear Dimensionality Reduction"
%   
%
%   Please ensure you have an up-to-date MATLAB version (> 2017a) as older
%   versions have a problems handling character/string arrays in certain 
%   cases which are extensively used in this script.
%
%   The script is segmented into four main categories:
%
%    -- Synthetic data evaluation: bench PM, MOSES, FD, RFD, & GROUSE using  
%                                  synthetic datasets
%    -- Real data evaluation: bench PM, MOSES, FD, RFD, & GROUSE using real 
%                             datasets
%    -- Speed tests: compare the execution speed of MOSES when compared 
%                    with PM, FD, RFD, & GROUSE
%    -- MOSES scaling tests: compare the performance of MOSES, in terms of
%                            error across different parameters of 
%                            block size (b), rank (r), and ambient dim. (n)
%
% Author: Andreas Grammenos (ag926@cl.cam.ac.uk)
%
% Last touched date: 30/12/2018
% 
% License: 
%  code: GPLv3, author: A. Grammenos 
%  paper: A. Eftekhari, R. Hauser, and A. Grammenos retain their respective 
%         copyrights (pre-print link: https://arxiv.org/abs/1806.01304)
%
%           

%% Initialisation

% clear/close everything
clc; clear; close all;

% enable for reproducibility, comment for (slightly) different 
% (~random) results
rng(200);

% declare global variables
global pflag
global datasetPath
global use_fast_moses_only
global use_offline_svds
global use_fdr
global use_blk_err
global pdf_print
global fig_print

global run_synthetic
global run_real
global run_speed_test
global run_moses_scaling

global run_full_scaling
global run_exp1
global run_exp2
global run_exp3

% experiments to run
run_synthetic = 1;      % run synthetic evaluation (set 0 to skip)
run_real = 0;           % run real data evaluation (set 0 to skip)
run_speed_test = 0;     % run the calc. speed tests (set 0 to skip)
run_moses_scaling = 0;  % run the scaling moses tests (set 0 to skip)

% global flags setup

% printing flags
pflag = 1;              % print resulting figures to ./graphs/
pdf_print = 0;          % print resulting figures as .pdf
fig_print = 1;          % print resulting figures as .fig

% execution configuration
use_fast_moses_only = 1;% speed up by using fast moses <-- USE IT :)
use_offline_svds = 0;   % drastically speed up execution by disabling 
                        % offline svds calculation WARNING THIS OPTION IS
                        % PAINFULLY SLOW. <- DEF. DISABLE IT :)
use_fdr = 0;            % use robust fd -- same as fd but on the recon.
                        % we normalise using a*Id; using the shifted
                        % subspace by a*Id does not work well in our case.
use_blk_err = 0;        % calc. errors per block not per column
                        % provides a DRASTIC improvement in speed but less
                        % granular error reporting. For GROUSE & FD is 100
                        % for PM and MOSES is equal to their respective 
                        % block sizes for each run. <- Prob. use it
                        
% moses scaling flags
run_full_scaling = 0;   % no need to run unless performing full error check
run_exp1 = 1;           % run experiment 1: fixed b, r, variable n
run_exp2 = 1;           % run experiment 2: fixed r, n, variable b
run_exp3 = 1;           % run experiment 3: fixed b, n, variable r

% setup the vars for the workspace
setup_vars();

%% Synthetic data

if run_synthetic == 1
  fprintf("\n ** Running synthetic data evaluation ** \n");

  % Define synthetic data problem parameters
  n = 200;        % Ambient dim
  r = 10;         % Algorithm aims to find rank-r truncation of input data
  T = 10*n;       % Scope of the algorithm (namely, max time)
  nSim = 10;      % number of simulations

  % NOTE: change alpha to see differences.
  alpha = [0.01, 0.1, 0.5, 1];
  
  % run the evaluation loop
  for i = 1:size(alpha, 2)
      synthetic_dataset_eval(n, T, r, alpha(i), nSim);
  end
  
  fprintf("\n ** Finished synthetic data evaluation ** \n");
else
  fprintf("\n ** Skipping synthetic data evaluation **\n");
end

%% Real data

% check if we run real data
if run_real ~= 1
  fprintf("\n ** Skipping real data evaluation **\n");
else
  fprintf("\n ** Running real data evaluation ** \n");

  % setup path datasets
  lightData = strcat(datasetPath, 'q8calibLight.dat');
  tempData = strcat(datasetPath, 'q8calibHumTemp.dat');
  voltData = strcat(datasetPath, 'q8calibVolt.dat');
  humidData = strcat(datasetPath, 'q8calibHumid.dat');

  % Light data
  light_r = 20; % target rank for light data
  real_dataset_eval(lightData, light_r, "Light");

  % Temperature data
  temp_r = 20;  % target rank for temperature data
  real_dataset_eval(tempData, temp_r, "Temperature");

  % Volt Data
  volt_r = 20;  % target rank for voltage data
  real_dataset_eval(voltData, volt_r, "Voltage");

  % Temperature data
  humid_r = 20; % target rank for humidity data
  real_dataset_eval(humidData, humid_r, "Humidity");
  
  fprintf("\n ** Finished real data evaluation ** \n");
end

%% Speed tests

if run_speed_test ~= 1
  fprintf("\n ** Skipping algorithm speed evaluation **\n");
else
  fprintf("\n ** Running algorithm speed evaluation **\n");
  
  % power law distribution params
  alpha = 1;  
  % no. of trials
  trials = 5;
  % different ambient dims
  n_arr = 200:200:1000;
  
  % After setting the parameters, run the speed tests
  
  r = 1;  % target rank
  fprintf("\n !! Testing thin-r recovery n >>>> r, with r=%d !!\n", r);
  speed_test(n_arr, r, alpha, trials)
  
  r = 10; % target rank
  fprintf("\n !! Testing avg-r recovery n >>> r, with r=%d !!\n", r);
  speed_test(n_arr, r, alpha, trials)
  
  r = 50; % target rank
  fprintf("\n !! Testing fat-r recovery n > r, with r=%d !!\n", r);
  speed_test(n_arr, r, alpha, trials)
  
  r = 100; % target rank
  fprintf("\n !! Testing super fat-r recovery n > r, with r=%d !!\n", r);
  speed_test(n_arr, r, alpha, trials)
  
  
  fprintf("\n ** Finished algorithm speed evaluation **\n");
end

%% MOSES Scaling tests

if run_moses_scaling ~= 1
  fprintf("\n ** Skipping MOSES scaling evaluation **\n");
else
  fprintf("\n ** Running MOSES scaling evaluation **\n");
  
  % Scaling test parameters
  
  n_arr = 200:200:1200; % ambient dimension array
  r_arr = 5:5:25;       % r-rank
  m_blk_mul = 1:1:15;   % block multiplier (we are bound by r)
  
  % Execute the scaling test
  moses_scaling(n_arr, r_arr, m_blk_mul);
  
  fprintf("\n ** Finished MOSES scaling evaluation **\n");
end

%% Comparison script end.