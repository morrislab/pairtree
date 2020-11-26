Pairtree
========
Pairtree infers the phylogeny underlying a cancer using genomic mutation data.
Pairtree is particularly suited to settings with multiple tissue samples from
each cancer, providing separate estimates of subclone frequency from each sample
that constrain the set of consistent phylogenies.  The Pairtree algorithm is
described in {insert link to paper}. The algorithm consists of two phases:

1. Compute pairwise relation tensor over all mutations (or clusters of
   mutations). This provides the probability over each of four possible
   evolutionary relationships between each pair of mutations `A` and `B`.

2. Use the pairwise relation tensor to sample possible phylogenies, assigning a
   likelihood to each. As each phylogeny is sampled, Pairtree computes the
   subclone frequency of each mutation (or cluster of mutations) within the
   tree, balancing the need to fit the observed mutation data while still
   obeying tree constraints.

The algorithm is described in [Pairtree: fast reconstruction of cancer
evolutionary history using pairwise mutation
relationships](https://www.biorxiv.org/content/10.1101/2020.11.06.372219v1)
(Wintersinger et al.).


Installing Pairtree
===================
1. Install dependencies. To ease installation, you may wish to use
   [Anaconda](https://www.anaconda.com/products/individual), which includes
   recent versions of Python 3, NumPy, and SciPy.
   [Miniconda](https://docs.conda.io/en/latest/miniconda.html) also works well.
   You need the following:

    * Python 3.6 or greater
    * NumPy
    * SciPy
    * [tqdm](https://github.com/tqdm/tqdm) (e.g., install via `pip3 install --user tqdm`)
    * [Numba](https://numba.pydata.org/) (e.g., install via `pip3 install --user numba`)
    * [colorlover](https://github.com/plotly/colorlover) (e.g., install via `pip3 install --user colorlover`)
    * C compiler (e.g., GCC)

   Pairtree has only been tested on Linux systems, but should work on any
   UNIX-like OS (including macOS).

2. Clone the Pairtree repository, then download and build the C code required
   to fit subclone frequencies to the tree. This algorithm was published in
   [Jia et al.](https://arxiv.org/abs/1811.01129), and uses the authors' implementation with minor modifications.

        git clone https://github.com/jwintersinger/pairtree
        cd pairtree/lib
        git clone https://github.com/jwintersinger/projectppm
        cd projectppm
        bash make.sh

3. Test your Pairtree installation.

        cd ../../example/
        mkdir results && cd results
        # Run Pairtree.
        ../../bin/pairtree --params ../example.params.json ../example.ssm example.results.npz
        # Plot results in an HTML file.
        ../../bin/plottree --runid example ../example.ssm ../example.params.json example.results.npz example.results.html
        # View the HTML file.
        firefox example.results.html


Interpreting Pairtree output
============================
(add note about how logs will be written in JSON format if stdout/stderr is directed to a file)

Input files
===========
(document the input files)


Tweaking Pairtree options
=========================
Tweaking some of Pairtree's options can improve the quality of your results or
reduce the method's runtime. The most relevant options are discussed below.
For a full listing of tweakable Pairtree options, including algorithm
hyperparameters, run `pairtree --help`.

Finding better trees by running multiple MCMC chains
----------------------------------------------------
Pairtree samples trees using MCMC chains. Like any optimization algorithm
working with a non-convex objective, Pairtree's chains can become stuck in
local optima. Running multiple MCMC chains helps avoid this, since each chain
starts from a different point in space, and takes a different path through that
space.

By default, Pairtree will run a separate MCMC chain for each CPU core (logical
or physical) present. As these chains are run in parallel, using multiple
chains does not increase runtime. For more details, please refer to the
"Running computatins in parallel to reduce runtime" section below.

Changing number of samples, burn-in, and thinning
-------------------------------------------------
Three options control the behaviour of each MCMC chain used to sample trees.

Changing number of samples, burn-in, and thinning
-------------------------------------------------
Three options control the behaviour of each MCMC chain used to sample trees.

* Number of MCMC samples: by changing the `--trees-per-chain` option, you can
  make each chain sample more or fewer trees. The more samples each chain
  takes, the more likely it is that those samples will be from a good
  approximation of the true posterior distribution over trees permitted by your
  data. By default, Pairtree will take 3000 total samples per chain.

* Burn-in: each chain will take some number of samples to reach high-density
  regions of tree space, such that those trees come from a good approximation
  of the true posterior. The `--burnin` paramemter controls what proportion of
  trees Pairtree will discard from the beginning of each chain, so as to avoid
  including poor-quality trees in your results. By default, Pairtree discards
  the first one-third of samples from each chain.

* Thinning: MCMC samples inherently exhibit auto-correlation, since each
  successive sample is based on a perturbation to the preceding one. We can
  reduce this effect by taking multiple samples for each one we record. You can
  change the `--thinned-frac` option to control this behaviour. By default,
  `--thinned-frac=1`, such that every sample taken is recorded. By changing
  this, for instance, to `--thinned-frac=0.1`, Pairtree will only record every
  tenth sample it takes. By recording only a subset of samples taken, the
  computational burden associated with processing the results (e.g., by
  computing summary statistics over the distribution of recorded trees) is
  reduced, alongside the storage burden of writing many samples to disk.

Running computations in parallel to reduce runtime
--------------------------------------------------
Pairtree can leverage multiple CPUs both when computing the pairwise relations
tensor and when sampling trees. Exploiting parallelism when computing the
tensor can drastically reduce Pairtree's runtime, since every pair can be
computed independently of every other. Conversely, exploiting parallelism when
sampling trees won't necessarily improve runtime, since trees are sampled using
MCMC, which is an inherently serial process, While no single MCMC chain can be
run in parallel, however, you can run multiple independent chains in parallel,
as discussed in the previous section. This will help Pairtree avoid local
optima and produce better-quality results.

By default, Pairtree will run in parallel on all present CPUs. (On
hyperthreaded systems, this number will usually be inferred to be twice the
number of physical CPU cores.) You can change this behaviour by specifying the
`--parallel=N` option, causing Pairtree to use `N` parallel processes instead.
When computing the pairwise tensor, Pairtree will compute all pairs in
parallel, such that computing the full tensor using `N` parallel processes
should only take `1/N` as much time as computing with a single process.  Later,
when sampling trees, Pairtree will by default run as many independent MCMC
chains as there are parallel processes. In this scenario, if you have `N` CPUs,
running only a single chain will be no faster than running `N` chains, since
each chain executes concurrently. However, in the `N`-chain case, you will
obtain `N` times as many samples, improving your result quality.

You can explicitly specify `--tree-chains=M` to run `M` chains instead of the
default, which is to run a separate chain for each parallel process. If `M` is
greater than the number of parallel processes `N`, the chains will execute
serially, increasing runtime. (Given the serial execution, you should expect
that tree sampling with `M > N` will take `ceiling(M / N)` times as long as

Fitting subclone frequencies to trees
-------------------------------------
Critical to the Pairtree algorithm is the step of fitting subclone frequencies
to trees. Beyond the user's interpretation of subclone frequencies for
downstream analysis of Pairtree's results, the subclone frequencies will also
affect tree structure -- the data likelihood for a tree is determined by how
well the tree-constrained subclone frequencies fit the observed VAF data. We
provide several different algorithms for computing these frequencies that
strike different balances between computational speed and result accuracy.

* `--phi-fitter=projection`: By default, Pairtree uses the "Efficient
  Projection onto the Perfect Phylogeny Model" algorithm developed in [Jia et
  al.](https://arxiv.org/abs/1811.01129) This uses a Gaussian approximation of
  the binomial likelihood we're trying to maximize.

* `--phi-fitter=rprop`: Uses the [rprop variant of gradient
  descent](https://ieeexplore.ieee.org/document/298623). Separate step sizes
  are maintained for each frequency scalar in each tissue sample. Those step
  sizes are increased when the direction of a step is consistent with the
  previous step taken, and reduced otherwise. Generally, this algorithm will
  produce more accurate results than `projection` (i.e., higher data
  likelihoods for the same tree), since it is directly optimizing the binomial
  objective rather than a Gaussian approximation, but at a higher computational
  cost. The extra computational burden becomes greater with more subclones
  (e.g., 30 or more), while being insignificant for small trees (e.g., with ten
  subclones or fewer).

* `--phi-fitter=proj_rprop`: First run the default `projection` algorithm, then refine its
  results using `rprop`. While this can produce better subclone frequency values
  at the cost of more computation, it has in testing demonstrated little
  benefit to accuracy.

* `--phi-fitter=graddesc`: uses gradient descent with an adaptive step size.
  (I.e., the step size is increased each time a step is accepted, and reduced
  each time a step is rejected.) Only a single global step size is used for all
  subclone frequencies across subclones and samples.
