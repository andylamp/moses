# Introducing MOSES

MOSES is a new take at tackling the problem of (M)emory-limited 
(O)nline (S)ubspace (Es)timation; in particular, MOSES can estimate 
the principal components of data and also reduce its dimension. In 
terms of its origins, MOSES slightly generalises the popular 
incremental Singular Vale Decomposition (iSVD) to work with 
*thin blocks* of data. This simple generalisation is in part what 
allows us to complement MOSES with a comprehensive statistical 
analysis that was not available for incremental SVD, despite its 
empirical success (more [here][5]). This generalisation also 
enables us to concretely interpret MOSES as an approximate solver 
for the underlying non-convex optimisation problem. We also find 
that MOSES shows state-of-the-art performance in our numerical 
experiments with both synthetic and real-world datasets while 
being orders of magnitude faster than its competitors when the 
ambient dimension (`n`) is large (>1000).

# Requirements

The code is generally self contained and all datasets are included or 
generated thus, in theory, just having `Matlab`  installed should be more than 
enough. It has to be noted though that due the recent `Matlab` changes on 
how it handles character and string arrays you should use a recent 
version of it -- the code was developed and tested in `Matlab` `2017b` and 
tested also on versions `2017a`, `2018a`; moreover, to address different 
OSes, care has been taken so that this code runs without any problems both 
on Windows-based machines as well as Unix-based ones.

# Streaming, memory limited, r-truncated SVD Method Comparison

In this instance we perform a comparison using both synthetic and real 
data against three methods which compute an approximate *memory-limited, 
streaming r-truncated SVD*. These methods are the following:

 * MOSES (https://arxiv.org/pdf/1806.01304.pdf)
 * Power Method (https://arxiv.org/pdf/1307.0032.pdf)
 * GROUSE (https://arxiv.org/pdf/1702.01005.pdf)

# Running the comparison

Running the comparison is simple -- just `cd` to the cloned `moses` directory 
within `Matlab` and run `comparison.m`. Running might take a while, if you want
to speed things up just try altering the setup parameters shown below:

```Matlab
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
use_offline_svds = 1;   % drastically speed up execution by disabling 
                        % offline svds calculation WARNING THIS OPTION IS
                        % PAINFULLY SLOW. <- DEf. DISABLE IT :)
use_blk_err = 0;        % calc. errors per block not per column
                        % provides a DRASTIC improvement in speed but less
                        % granular error reporting. For GROUSE it is 100
                        % for PM and MOSES is equal to their respective 
                        % block sizes for each run. <- Prob. use it
                        
% moses scaling flags
run_full_scaling = 0;   % no need to run unless performing full error check
run_exp1 = 1;           % run experiment 1: fixed b, r, variable n
run_exp2 = 1;           % run experiment 2: fixed r, n, variable b
run_exp3 = 1;           % run experiment 3: fixed b, n, variable r
```

Please note that you can tweak the relevant section values if you want 
to run slightly different experiments but if you want to reproduce the 
results in the paper please leave these values as-is.

# Synthetic datasets

The synthetic dataset is measured using random vectors drawn from a power
law distribution with the following alpha values in this instance: `0.01`, 
`0.1`, `0.5`, and `1` while lambda always set to `1`. Practically speaking
this is eloquently materialised by using the following segment:

```Matlab
% generate the singular spectrum
dd = lambda*(1:n).^(-alpha);
% generate Sigma
Sigma = diag(dd);
% random initialization of S basis
S = orth(randn(n));
% given S and Sigma generate the dataset (Y)
Y = (S * Sigma * randn(n, T))/sqrt(T-1);
```

# Real datasets

The real datasets are the the ones supplied with [this][3] paper 
retrieved from [here][1] and they are the following:

 * Light Data (48x7712)
 * Humidity Data (48x7712)
 * Volt Data (46x7712)
 * Temperature Data (56x7712)

# Error metrics

To compare MOSES against Power Method and GROUSE we employ the following 
two metrics:

 * The Frobenius norm of `Yr` vs `Y` columns seen so far normalised using 
    their respective arrival time.
 * The final Frobenius norm of `Yr` vs `Y` normalised by the final `T`.

For speeding up the computation there are two modes that calculate the 
over-time errors. The first one is *per-block* while the second 
is *per-column* the latter being significantly slower than the former. 
For completeness, in our paper we used the *per-column* error 
calculation throughout but similar error metrics can be achieved by 
using the faster *per-block* errors.

Toggling this particular mode is done by setting `use_blk_err`
in `comparison.m` to either `1` (on) / `0` (off).

## Normalised Frobenius norm over time normalised with current T

The error metrics are calculated using the Frobenius norm for the
matrix columns seen so far normalised by the current time. The full 
formula to find the error at column `k` would be:

```Matlab
ErrFro(k) = sum(sum((Y(:, 1:t)-YrHat_c).^2, 1))/t;
```

Where `YrHat_c` is:

```Matlab
SrHatTemp = SrHat(:, 1:r); % r-truncation of the SVD
% SrHat in this instance is the previous block subspace estimation
YrHat_c = (SrHat*SrHat')*Y(:, 1:k*B); 
```

## Final Yr vs Y Frobenius norm over T

The other metric is the Frobenius norm difference of the three `Yr` 
approximation methods vs. the original values `Y` 
normalised using the total number of columns seen (`T`).

The formula is:

```Matlab
froT = n * Sum_{1}_{n} [ (Y_i - Yr_i)^2 ] / T
```

This metric shows us another view of how different these resulting 
matrices are when compared to the original ones.

## Plots

A number of plots are generated while running the comparison and for
convenience they are printed into a generated directory under the 
`graph` directory. Each directory is named using the current 
timestamp upon creation as its name and the timestamp format 
following the [ISO-8601][4] standard.

Additionally, the printing function is flexible enough to able to export 
in three commonly used formats concurrently -- namely `png`, `pdf`, and 
`fig` for easier processing. Of course, by toggling the appropriate flags 
printing to `pdf` and `fig` can be disabled thus saving space. For 
brevity these are the following:

```MatLab
% printing flags
pflag = 1;              % print resulting figures to ./graphs/
pdf_print = 0;          % print resulting figures as .pdf
fig_print = 1;          % print resulting figures as .fig
```

### Synthetic data

A number of different plots are generated for the synthetic data, these
are explained below: 

 * The first one is a `subplot` containing two plots; the first one shows 
   the scree plot of the singular values while the second one shows the 
   *mean* Frobenius normalised error for the three methods over time for 
   the total number of simulations run for that particular *alpha* value.

 * The second plot is the second `subplot` from above without the scree
   plot.

 * The third plot shows the averaged error over time for the number 
   of trials performed but just for Power Method and MOSES; this is done in 
   order to have a better comparison as GROUSE results sometimes 
   skewed the graph.

 * The fourth plot is again a `subplot` which the first one shows the 
   resulting final error per iteration for Offline SVD, MOSES, 
   Power Method, and GROUSE; while the second one shows the error for 
   just MOSES and Power Method as they are closer.

 * The fifth, and final, plot for the synthetic data shows the final
   error per iteration for only the Power Method and MOSES.

### Real data

We only plot the Frobenius normalised error in a streaming fashion as 
the columns are revealed for `Yr` vs `Y` to the three methods 
over time for each of the four datasets described previously;
for readability, all figures are generated in a logarithmic scale 
using `semilogy`. Before processing the datasets it is important to note
that the data are centered using the following:

```Matlab
% update the column number
cols = size(Y, 2);
% center the dataset by using Y = Y - ((Y*(vec_ones*vec_ones'))./cols)
vec_ones = ones(cols, 1);
Y = Y - ((Y*(vec_ones*vec_ones'))./cols);
```

### Speed tests

For speed tests we tests how the algorithms perform as the ambient 
dimension (`n`) and rank (`r`) change. For our experiments we used 
three rank (`r`) values: `r = [1, 10, 50]` and an ambient dimension
range: `n = 200:200:1200`. In this instance three plots are 
generated for the tests, which are the following:

 * Speed for all three method when using:`r=1` and `n = 200:200:1200`
 * Speed for all three method when using:`r=10` and `n = 200:200:1200`
 * Speed for all three method when using:`r=50` and `n = 200:200:1200`

### MOSES scaling tests

To evaluate how MOSES scales with respect to block size (`b`), 
rank (`r`), and ambient dimension (`n`) we performed three tests
in which we fixed two out of the three parameters and plotted
the average error out of `10` trials. The three resulting 
plots are the following:

 * Scaling experiment 1: fixed `r=15`, `b=2r`, and `n = 200:200:1200`
 * Scaling experiment 2: fixed `r=15`, `n=1200`, and `b= r:1:15r`
 * Scaling experiment 3: fixed  `n=1200`, `b=2r`, and `r= 5:5:25`

# Code Organisation

The code is organised mainly in the following files:

 * `comparison.m`: The main starting point of the experiments, initial parameters are defined there.
 * `grouse.m`: Original `GROUSE` algorithm code as provided from the paper
 * `mitliag_pm.m`: Implementation of Mitliagkas Power Method for Streaming PCA
 * `moses_fast.m`: A more efficient implementation of MOSES
 * `moses_scaling.m`: function that is responsible for running the scaling experiments for MOSES
 * `moses_simple.m`: A simple implementation of MOSES
 * `my_grouse.m`: Wrapper to run `grouse.m` which sets execution parameters (as seen [here][2])
 * `my_toc.m`: function that processes the `toc` with better formatting
 * `online_svds_real.m`: function which compares the three methods against a real dataset
 * `online_svds_synthetic.m`: function which compares the three methods against a synthetic dataset
 * `print_fig.m`: function that, if enabled, outputs the graphs generated for current run as images or `pdf`
 * `README.md`: This file, a brief README file.
 * `real_dataset_eval.m`: function which is a wrapper for `online_svds_real` that processes its results and draws the plots
 * `setup_vars.m`: function which sets the correct variables based on OS
 * `synthetic_data_gen.m`: function which generates a matrix with random vectors from a power law distribution
 * `synthetic_dataset_eval.m`: function which is a wrapper for `online_svds_synethic` that processes its results and draws the plots

 # License

This code is licensed under the terms and conditions of GPLv3 unless 
otherwise stated. The actual paper is governed by a separate license and 
the paper authors retain their respective copyrights.

# Acknowledgement

If you find our paper useful or use this code, please consider citing our work 
as such:

```
@misc{1806.01304,
Author = {Armin Eftekhari and Raphael A. Hauser and Andreas Grammenos},
Title = {MOSES: A Streaming Algorithm for Linear Dimensionality Reduction},
Year = {2018},
Eprint = {arXiv:1806.01304},
}
```

# Disclaimer

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS “AS IS” 
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR 
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES 
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; 
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND 
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS 
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

[1]: http://www.cs.cmu.edu/afs/cs/project/spirit-1/www/
[2]: http://web.eecs.umich.edu/~girasole/grouse/
[3]: http://www.cs.albany.edu/~jhh/courses/readings/desphande.vldb04.model.pdf
[4]: https://en.wikipedia.org/wiki/ISO_8601
[5]: https://arxiv.org/pdf/1806.01304.pdf


