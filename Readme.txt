Poisson embeddings
------------------

This implements exponential family embeddings for shopping data, as described in

  M. Rudolph, F. J. R. Ruiz, S. Mandt, and D. M. Blei. Exponential Family Embeddings.
  Advances in Neural Information Processing Systems. Barcelona (Spain), December 2016.
  [http://papers.nips.cc/paper/6571-exponential-family-embeddings.pdf]

USAGE:
	pemb [options]

OPTIONS:

-dir <string>          path to directory with all the data files (see format below)
                       [default: .]

-outdir <string>       path to directory where output files will be written
                       [default: .]
-userVec <int>         use per-user vectors (if 0, no per-user vectors are
                       included; if 1, the per-user vectors are added to the
                       context; if 2, the per-user vectors are multiplicative;
                       if 3, the per-user vectors interact with alpha)
                       [default: 0] 

-K <int>               number of latent factors
                       [default: 50]
-itemIntercept         use per-item intercepts
                       [default: disabled]
-userVec <int>         use per-user vectors (if 0, no per-user vectors are
                       included; if 1, the per-user vectors are added to the
                       context; if 2, the per-user vectors are multiplicative;
                       if 3, the per-user vectors interact with alpha)
                       [default: 0] 
-avgContext <int>      define context averaging (if 0, no averaging is
                       performed; if 1, the context is averaged by the sum of
                       ratings/units; if 2, the context is averaged by the number
                       of unique movies/items)
                       [default: 1]
-binarizeContext       the observables in the context are treated as binary,
                       instead of using the actual counts
                       [default: disabled] 
-additiveMean          specify that the Poisson mean is directly the inner
                       product, without exponentiation
                       [default: disabled]
-gaussian              use Gaussian likelihood instead of Poisson
                       [default: disabled]
-bernoulli             use Bernoulli likelihood instead of Poisson
                       [default: disabled]

-batchsize <int>       number of baskets to form minibatches (if <=0, uses all
                       baskets in the dataset)
                       [default: -1]
-negsamples <int>      number of negative samples (if <=0, all negative samples
                       are included in the objective function)
                       [default: 10]
-nsFreq <int>          choose the distribution for negative sampling (-1: uniform;
                       0: unigram; 1: unigram raised to the power of 3/4)
                       [default: -1]
-zeroFactor <double>   factor (0<zeroFactor<=1) to downweight the zeroes (if <=0,
                       zeroFactor is set to 1)
                       [default: 0.1]
-max-iterations <int>  maximum number of iterations
                       [default: 1500]

-keepOnly <int>        only consider the N most frequent items; ignore the rest
                       (in case of ties, this can keep more than N items). If <=0;
                       do not ignore any item
                       [default: -1]
-keepAbove <int>       only consider the items that appear in at least N baskets;
                       ignore the rest (if <=0; do not ignore any item)
                       [default: -1]

-noVal                   ignore the validation file (if present) and stop the
                         simulation only after max-iterations
                         [default: disabled ---the validation set is used]
-valTolerance <double>   relative tolerance of the validation log-likelihood to
                         declare convergence
                         [default: 0.000001]
-valConsecutive <double> number of consecutive iterations that validation log-
                         likelihood is allowed to decrease before declaring
                         convergence
                         [default: 5]
-noTest                  ignore the test file (if present)
                         [default: disabled ---the test set is used and
                         log-likelihood is reported to a file]

-adagrad               use Adagrad for stochastic optimization (default)
                       [default: enabled]
-rmsprop               use RMSProp instead of Adagrad
                       [default: disabled]
-eta <double>          stepsize for the optimization algorithm
                       [default: 0.1]
-gamma <double>        forgetting factor for the optimization (only if RMSProp)
                       [default: 0.9]

-s2rho <double>        variance of the prior over the embedding vectors
                       [default: 1.0]
-s2alpha <double>      variance of the prior over the context vectors
                       [default: 1.0]
-s2theta <double>      variance of the prior over the user vectors
                       [default: 1.0]
-s2all <double>        set all the prior variances
-s2noise <double>      variance of the Gaussian likelihood (only if Gaussian
                       model)
                       [default: 1.0]
-noRegulariz           deactivate regularization
                       [dafault: disabled ---regularization is performed]

-stdIni <double>       standard deviation used only for random initialization of
                       the parameters
                       [default: 0.1/sqrt(K) for the non-additive model;
                        0.1 for the additive model]
-iniPath <string>      path to folder that contains the files to initialize
                       all the parameters
                       [default: disabled ---random initialization]
-iniThetaVal <double>  initialize all user vectors to this value (overwrites
                       '-iniPath' for the user vectors).
                       [default: disabled ---random initialization]

-seed <int>            set random seed
                       [default: 0]
-rfreq <int>           frequency (#iterations) for evaluating convergence and
                       computing predictive performance
                       [default: 10]
-saveCycle <int>       frequency (#iterations) to save output files
                       [default: 100]
-label <string>        label to be appended to the name of the output folder
                       [default: none]


Compilation
-----------

You need to compile the code using C++11 and linking to the GSL libraries. 
Here are some examples:

  g++ -std=c++11 -Wall -o pemb pemb.cpp `gsl-config --cflags --libs`
  g++ -Wall -I/usr/local/include -lgsl -o pemb pemb.cpp
  g++ -std=c++11 -Wall -lgsl -lgslcblas -o pemb pemb.cpp


Input
-----

The model needs the following input files:

train.tsv       contains the training data. It has three columns, in
                tab-separated format, without headers. First column:
                user id (a non-negative integer). Second column: item id
                (a non-negative integer). Third column: session id (a
                non-negative integer). Fourth column: number of units
                (a non-negative integer). Thus, the format of each line of
                train.tsv is as follows (tab-separated):
                    user_id   item_id   session_id   units
                The session id allows to specify different trips for each
                user. Session ids for two different users are independent
                (i.e., user 1 may make two trips, with session ids 1 and 2,
                while user 2 may make only two trips, with session ids 2 and
                3. The program will find that each user makes two trips,
                regardless of the specific numbers used for the session ids)
                IMPORTANT: Do not include the zeros in any of the input files
                (only those transactions with units>=1). Also: The code assumes
                that data is sparse. If this is not the case, the program will
                run significantly slower.
test.tsv        contains the test data, in the same format as train.tsv. The
                test data should contain new sessions, independent of the
                sessions in train.tsv. In other words, to compute the performance
                on the test set, the context will be defined by the information
                present in test.tsv *only*
validation.tsv  contains the validation set, used to assess convergence, in the
                same format as train.tsv. In contrast to test.tsv, the context
                will be defined by the information present in train.tsv *only*.
                In other words, for each line of validation.tsv, the program will
                look for all the lines with the same (user_id,session_id) pair
                in train.tsv in order to obtain the context

Output
------

+param_rho.txt:   The embedding vectors
+param_alpha.txt: The context vectors.
                  Both param_rho.txt and param_alpha.txt are tab-separated
                  TSV files (with no header), in which each line contains:
          line_number    item_id    value_1    [value_2    [value_3    [...]]]

+param_theta.txt: If '-userVec N' is specified with N>0, this file will also
                  be created, containing the per-user vectors. Its format is:
          line_number    user_id    value_1    [value_2    [value_3    [...]]]

NOTE: Intermediate files (saved every saveCycle iterations) will also be created
      with the same format. The file names will contain the iteration number.

+test.txt: Performance on the test set.
+validation.txt: Performance on the validation set.
                 These are tab-separated TSV files (without header) containing
                 the following information:
          iteration    duration    llh1    llh2    llh3    n    ny
                 Here, 'iteration' is the iteration number (reported only
                 every '-rfreq' iterations), 'duration' is the total wall-clock
                 duration of the program (in seconds) since it was started,
                 'llh1' is the Poisson log-likelihood, 'llh2' is the log-likelihood
                 after normalizing the Poisson mean for each possible item, 'llh3'
                 is the same as 'llh2' but it is "binary" in the sense that it ignores
                 the number of units of each item (or the ratings of movies) and only
                 considers their presence, 'n' is the number of items used to
                 compute these numbers, and 'ny' is the sum of the units (or ratings)
                 for these items: 'n' and 'ny' should remain constant across all
                 iterations because they depend on the information in test/validation
                 only

+log.txt: a log file containing a list of the parameters used in the experiment
+max.txt: a file containing (in order) the final number of iterations, the duration
          of the experiment in seconds, the validation log-likelihood, and the reason
          why the program stopped (0=validation llh has converged; 1=validation llh
          is decreasing for several consecutive iterations; 2=reached maximum number
          of iterations)
+avgNorm.txt: Contains the average norm of the latent vectors at each iteration. In
              order: iteration number, norm of rho, norm of alpha, norm of theta (the
              last is only present if '-userVec')
+obj_function.txt: The value of the objective function (stochastically estimated).
                   In order: iteration number, objective function
+telapsed.txt: Wall-clock time elapsed per iteration. In order: iteration number,
               elapsed time per iteration (in seconds)
