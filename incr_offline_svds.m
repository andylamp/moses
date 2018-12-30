function [OffT, OffFroError, Yr, t] = incr_offline_svds(Y, r, floorMul, no_err)
%INCR_OFFLINE_SVDS Incremental offline svds of rank r, for comparison only
%don't use... very(, very!) inefficient.
%
% Author: Andreas Grammenos (ag926@cl.cam.ac.uk)
%
% Last touched date: 30/12/2018
% 
% License: GPLv3
%
  fprintf('\n ** Running (the painfully slow) Offline SVDS...\n');
  
  % scope in global variables
  global use_blk_err
  
  % default block size
  blk_size = 100;
  cnt = 1;

  % check if we have a multiplier
  if nargin < 3
    floorMul = 2;
  end
  
  % check if we have a no error flag
  if nargin < 4
    no_err = 0;
  end

  % get the number of rows
  numc = size(Y, 2);
  % set up YrHat_c
  YrHat_c = [];
  % start the clock
  ts = tic;
  % set this to be equal to the number of columns
  if no_err == 0
    % check if we use a block error
    if use_blk_err == 1
      OffT = floor(numc/blk_size);
      OffFroError = zeros(1, OffT);
    else
      OffT = 1:numc;
      OffFroError = zeros(1, numc);
    end
  else
    OffT = 0;
    OffFroError = 0;
  end

  % run the svds for Y -- BEWARE, VERY SLOW
  for cur_t = 1:numc
    y_k = Y(:, 1:cur_t);
    % run the svds
    [Soff, ~, ~, ~] = svds(y_k, floorMul*r);
    % reduce the components to the r-truncation
    if cur_t > r
      Soff = Soff(:, 1:r);
    end

    % Incremental offline svds errors start

    if no_err == 0
      if use_blk_err == 1
        % block error, less granular but much, much faster
        if mod(blk_size, cur_t) == 0
          % block error, faster than full but still painfully slow
          YrHat_c = (Soff*Soff')*Y(:, 1:cur_t);
          temp = sum(sum((y_k-YrHat_c).^2, 1));
          OffFroError(cnt) = temp/cur_t;
          cnt = cnt + 1;
        end
      else
        % full error, even more painfully slow
        YrHat_c = (Soff*Soff')*Y(:, 1:cur_t);
        temp = sum(sum((y_k-YrHat_c).^2, 1));
        OffFroError(cur_t) = temp/cur_t;
      end
    end

    % Incremental offline svds errors end

  end

  % update the final estimate
  Yr = YrHat_c;
  % report (the sadly insane) execution time
  t = my_toc(ts);
  % spacing
  fprintf("\n");
end

