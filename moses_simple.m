function [T, ErrFro, U, Yr, t] = moses_simple(Y, r, blk_size, floor_mul, no_err)
%% This is a simple implementation of MOSES (https://arxiv.org/pdf/1806.01304.pdf) 
%
% Author: Andreas Grammenos (ag926@cl.cam.ac.uk)
%
% Last touched date: 30/12/2018
% 
% License: GPLv3
%      
  fprintf('\n ** Running MOSES Simple...\n');
  
  % scope in the global variable
  global use_blk_err
  
  % get Y details
  [dim, Ti] = size(Y); 
  
  % check we calculate the error (disabled for speed runs)
  if nargin < 5
    no_err = 0;
  end

  % check if we have a floor multiplier as an argument
  if nargin < 4
    floor_mul = 2;
  end

  % check if n < r or n == 1
  if dim == 1 || dim < r
    error(" ** ERR: Ambient dimension must be > 1 and r < n **");
  end
  
  % moses configuration
  
  % set the block, depending on argument
  if nargin < 3
    b = 2*r;
  else
    if blk_size < 2*r
      fprintf("\n !! WARN: Block size must be at least 2*n !!\n");
      b = 2*r;
    else
      b = blk_size;
    end
  end
  
  % check if Ti < b, in which case we cannot run it
  if Ti < b
    error("\n Block size must be lower than the number of columns");
  end
  
  K = floor(Ti/b);              % Number of blocks
  cnt = 1;                      % index for error block align 
  
  % preallocate based on no error run and block error
  if no_err == 0
    if use_blk_err == 1
      T = nan(1, K);                % T steps for error log
      ErrFro = nan(1, K);           % Fro normalised error with T
    else
      T = (b+1):Ti;                 % T steps for error log
      ErrFro = nan(1, size(T, 2));  % Fro normalised error with T
    end
  else
    T = 0;
    ErrFro = 0;
  end
  % generated YrHat
  YrHat = [];
  % output the block number
  fprintf([' ** Total number of blocks (k): %d ', ...
    'with Block size of: %d\n'], K, b);

  % start timing
  ts = tic;
  for k = 1:K
    Yk = Y(:,(k-1)*b+1:k*b);
    min_t = ((k-1)*b)+1;
    max_t = k*b;
    Yk_cur = [YrHat, Yk];
    [SrHat, ~, ~] = svds(Yk_cur, floor(floor_mul * r));
    SrHatTemp = SrHat(:, 1:r);

    % MOSES errors start
    if no_err == 0
      if k > 1
        if use_blk_err == 1
          % Now calculate the Block normalised errors
          YrHat_c = (SrHatTemp*SrHatTemp')*Y(:, 1:max_t);
          % Frobenius norm incremental error, per block located at 
          % k*b normalised with current T.
          temp = sum(sum((Y(:, 1:max_t)-YrHat_c).^2, 1));
          ErrFro(cnt) = temp/max_t;
          T(cnt) = max_t;
          cnt = cnt + 1;
        else
          % Now calculate the Block normalised errors
          YrHat_c = (SrHatTemp_old*SrHatTemp_old')*Y(:, 1:max_t);
          % Frobenius norm incremental error, per column in block located  
          % at (k-1)b:kb normalised with current T.
          temp = cumsum(sum((Y(:, 1:max_t)-YrHat_c).^2, 1));
          ErrFro(cnt:cnt+b-1) = temp(min_t:max_t)./(min_t:max_t);
          cnt = cnt + b;
        end
      end
    end

    % MOSES errors end

    SrHatTemp_old = SrHatTemp;
    
    % update Yr estimation
    if k > 1
      if no_err == 0
        YrHat = YrHat_c;
      else
        YrHat = (SrHatTemp*SrHatTemp')*Y(:, 1:k*b);
      end
    end
  end
  U = SrHatTemp;   % current U estimation
  Yr = YrHat;      % current YrHat estimation
  t = my_toc(ts);  % current trial execution time
end