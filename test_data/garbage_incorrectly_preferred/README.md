This is an instance where the `garbage` model is incorrectly preferred over `diff_branches` and `B_A`, both of which should be equally likely. To recreate the problem, do the following:

1. Apply `spawn_repl.patch` to `bin/removegarbage`.

2. Run `NUMBA_DISABLE_JIT=1 python3 ~/work/pairtree/bin/removegarbage --pairwise-results pairs.npz --prior 0.01 --max-garb-prob 0.5 mats08.ssm mats08.params.json out2.json`.

3. When the REPL spawns, run `import lh; pairwise._examine('s0', 's1', variants, logprior={'garbage': -np.inf}, _calc_lh=lh.calc_lh_quad)[0][14]`.

4. Observe:
  1. The first output is a matrix. Each row is a sample. The columns are put into three groups and separated by `NaN`.
    a. The first batch of columns is `(var_reads, total_reads, vaf, omega_v, vaf/omega_v)` for the first variant.
    b. The second batch is the same tuple for the second variant.
    c. The third batch are the log-evidence values for the models in order `garbage=0, cocluster=1, A_B=2, B_A=3, diff_branches=4`
  2. The second output is the log-evidences summed across samples. You can see `garbage` is highest (`-468`), followed by `diff_branches` (`-478`).
  3. The third output is the posterior computed from the log-evidences, which here is using a uniform prior across the five models

5. Notes:
  * If we scan down the last batch of columns, we see that `diff_branches` is the best model in `21` of `29` samples
  * In the remaining `8` samples, we prefer `garbage` in 2 of them and `cocuster` in 6 of them.
  * In the two samples where `diff_branches` loses to `garbage`, it loses by `29 nats` total, which is enough to drive the `10 nat` difference we see when summing across samples. Because of how we defined the pairwise model, `garbage` will never be less than half as likely as any of the other models. This means, at best, when `diff_branches` wins over garbage, it will accumulate an advantage of `log(2) = 0.69` nats in each sample. Given that it wins in 21 samples, this gives it at most an edge of `21 * log(2) = 14.6 nats` from those samples

6. Source of bug and possible solutions:
  * This is specifically an error relating to 1D quadrature via QUADPACK. Integrating via 1D Monte Carlo or a more naive 2D Monte Carlo produces a result that is usually more-or-less correct, but very high variance, such that the correct model doesn't always win.
  * This only occurs at high read depths. Multiplying the var, ref, and total read counts by 1/100 or even 1/10 results in the right model winning for 1D QUADPACK.
  * Intuitively, the failure at high depth makes sense. If the mode of the integrand spans a tiny part of the unit interval, then quadrature can miss it entirely. By contrast, the garbage integral doesn't use quadrature -- it computes two separate integrals over the binomial, one for each variant in the pair, since there's no dependence between `phi1` and `phi2` for it. Instead, it can just integrate over the binomial using the incomplete beta.
  * We know approximately what the integrand looks like when we plot it -- we know where all the mass is going to be. Ideally, we want to be able to tell our quadrature routine, "You need to focus on integrating over a bunch of tiny intervals concentrated around the MLE of the integrand. Don't worry about stuff far away from it." Practically, I don't know how to do that.
  * Seemingly, the problem occurs when most of the mass is concentrated right at the edges of the definite integral's limits of integration (e.g., when most of the mass occurs around phi=0 or phi=1). This makes sense -- maybe QUADPACK can't "find" this mass. To overcome this, we can either scale down the total read count so that the mass becomes more spread out, or we can try to compute the integral analytically for the special case when `V=0` or `V=T`. (Is it actually at `V=T`? Or is it when `V` is in the range `[omega*T - epsilon, T]`?) Quaid made some progress in figuring out how to do this, but it only works in the `A_B` and `B_A` cases (not `diff_branches`), and it requires `omega_A = omega_B`.
