\name{NEWS}
\title{spaMM News}
\encoding{UTF-8}
% cf ?news for the format and usage of this file
%
\section{CHANGES UP TO VERSION 3.6.0}{ 
\subsection{NEW FEATURES}{
\itemize{
\item (major) New function fitmv() to fit multivariate-response models, extending most of the features of fitme().
\item LRT() (aka anova.HLfit()) can now perform a bootstrap test for models differing (also) by their random effects. 
\item fitme(), and the new fitmv(), accept partially-fixed random coefficients parameters, as in 'fixed=list(ranCoefs=list("1"=c(NA,0.5,NA)))'
\item get_inits_from_fit() now gets also the estimates of the shape parameters of the negbin and COMPoisson families. 
\item New function map_ranef() to plot a map of a spatial random effect (using only graphic functions from base R, but other functions can use the returned grid of values).
\item Improved control of names of user-provided fixed lambda values (new error message may result).
\item poly() in random-coefficient terms is now handled.
\item New extractor residVar() for residual variances, with features different from get_residVar(). 
}
}
\subsection{USER-LEVEL CHANGES}{
\itemize{
\item Most GLMM fits using sparse-precision algorithm should be substantially faster (near 50\% for long fits).
\item Better control of progress bars as the global 'barstyle' option is no longer used in some cases where a local argument should be used.
\item The 'uniqueGeo' argument of HLCor() has long been obsolete and is now deprecated. 
\item VarCorr.HLfit() gets a new argument 'add_residVars' to make information about residual variance optional (also fixing potential bugs).
\item residuals() did not always return a vector.
}
}
\subsection{FIXED BUGS}{
\itemize{
\item confint() on PQL[/L] fits had several issues, in particular easily producing spurious messages suggesting that the fit was inaccurate.  
\item One could not use a Matern(RHS|.) (for non-unit RHS) together with a distMatrix for this Matern term.
\item Improve a backward incompatibility issue in the (internal) post-fit function determining which fitting function was used to produce a fit.
\item The recently introduced extractor get_RLRsim_args() did not correctly handle marginal LR tests (as opposed to REML ones).
\item the global 'barstyle' option was always directly used, so that the 'showpbar' argument did not control the progress bar in dopar().
\item plot.HLfit() reported incorrect standardized deviance residuals for binomial, non-binary response.
\item Several post-fit functions such as predict() could report incorrect results for 'IMRF' (and perhaps 'adjacency') models since version 3.4.1.
\item The (hardly ever tested) 'sqrt' family link was not fully implemented for families that allows it. 
\item summary() printed an incorrect number of groups for random-coefficients terms (long known and not bothered about until now). 
\item Documentation for the 'distMatrix' argument of HLCor() and 'higher' fitting functions was incomplete.
\item simulate(., re.form=NULL) was broken.
\item random-coefficient parameters fitted by "outer optimization" where not counted in the degrees of freedom for AIC.
\item Fitting some LMMs with fixed 'lambda' could return incorrect results. 
\item Typo in printed summary for rand.family=Gamma (wrong 'scale' value).
}
}
\subsection{ALL THE ABOVE for the 'multi' family}{
\itemize{
\item It was given a new scrutiny, leading to a display bug fix, a revised fitted() output, improved summary() information, a new how.HLfitlist() extractor, and revised documentation. But fitmv() might ultimately be more useful to fit multinomial data. % ___FIXME___
}
}
}
\section{CHANGES UP TO VERSION 3.5.0}{ 
\subsection{NEW FEATURES}{
\itemize{
\item 'Nested' Matern() or Cauchy() random effects are now fitted, using a syntax of the form Matern(1|x+y \%in\% grouping_factor).
\item 'barstyle' option controls more progress bars than previously.
\item New extractor get_RLRsim_args() extends previousy existing get_RLRTSim_args().
\item get_ranPars() handles new 'which' values for extraction of variance parameters.
}
}
\subsection{USER-LEVEL CHANGES}{
\itemize{
\item More judicious and limited production of warnings in some prediction variance computations.
}
}
\subsection{FIXED BUGS}{
\itemize{
\item An 'invalid type/length (S4/0) in vector allocation' bug occurring for recent spaMM versions on 3.6-ish versions of R.
\item Glitch in VarCorr() for conditional autoregressive fits (a line is now filled with NA's rather than with inappropriate values).
}
}
}
\section{CHANGES UP TO VERSION 3.4.0}{ 
\subsection{NEW FEATURES}{
\itemize{
\item \code{IMRF()} terms now also handle objects created by INLA::inla.spde2.matern() with non-integer alpha values, and therefore a broader range of smoothness values (PS: doc corrected wrt this only in version 3.6.0).
\item objects created by INLA::inla.spde2.pcmatern() are also handled.
\item 'sparse_precision' fitting method can be notably faster than before for large data sets. By default, this will affect fits of models with adjacency or IMRF terms. 
\item Fitting family=Gamma(log) mixed models will be faster than before, twofold or possibly more for large data sets.
\item confint.HLfit() now returns a list of confidence intervals when \code{parm} defines a vector of parameters. In the bootstrap case, the same bootstrap simulations can thus be used for several parameters.
\item get_predCov_var_fix() can now handle 'multIMRF' random effects correctly (for the happy few who would combine these two relatively cryptic functionalities).
\item \code{formula_env} specification to control a formula environment: see help("formula_env").
}
}
\subsection{USER-LEVEL CHANGES}{
\itemize{
\item The structure of the confint.HLfit() result is modified as this function no longer hacks the boot::boot.ci() results for printing purposes.
}
}
\subsection{FIXED BUGS}{
\itemize{
\item The 'distMatrix' argument was no longer effective. 
\item Prediction variance computations could fail for 'IMRF' random effects with 'newdata'.
\item confint() could fail when the object had been fitted by a "non-vanilla" call of the fitting function (e.g., 'spaMM::fitme' instead of 'fitme').
\item confint() by sparse-precision method could return incorrect results.
\item is_separated.formula() was not returning results properly for binomial non-binary response.
\item A case where LevenbergM for PQL could underfit was found.
\item confint(.,parm=<fixed-effect coefficient>) failed (stop) on sparse-precision fits with a single fixed-effect coefficient (typically "(Intercept)").
}
}
\subsection{DEPENDENCIES}{
\itemize{
\item spaMM is now formally declared to enhance 'INLA'.
}
}
}
\section{CHANGES UP TO VERSION 3.3.0}{ 
\subsection{NEW FEATURES}{
\itemize{
\item confint.HLfit() can now compute bootstrap confidence intervals by acting as a wrapper for boot::boot.ci(), and is then not limited to fixed-effect coefficients.
\item New hatvalues() extractor.
\item Faster and more memory-efficient step in Levenberg-Marquardt-like algorithm for "big" random effects fitted by sparse-correlation algorithms.
\item Faster and more memory-efficient step in prediction code for "big" random effects.
\item New option 'barstyle' to control globally the style of progress bars. New argument 'showpbar' for function dopar() and spaMM_boot().
\item Informative warning when some variables are not found in parallel context. New argument 'pretest_fn' for function dopar().
\item Better criteria for selection of sparse-precision algorithm, though far from perfect.
\item VarCorr() extractor should now handle all spaMM fit results.
\item Added 'seed' argument to spaMM_boot() and fixedLRT().
}
}
\subsection{DEPENDENCIES}{
\itemize{
\item 'ROI.plugin.glpk' moved from Imports to Suggests, removing a "hard" dependency to the external library glpk.
\item 'boot' in Imports.
}
}
\subsection{USER-LEVEL CHANGES}{
\itemize{
\item Added formal argument 'intervals' for get_intervals(). This does not affect any functionality, but explicitly shows that different intervals can be computed by changing this argument.  
\item All mixed-effect fitting functions now have a formal argument 'method' equivalent to 'HLmethod'. 'HLmethod' can still be used in a backward-compatible way.
}
}
\subsection{FIXED BUGS}{
\itemize{
\item simulate(<Gamma-response fit>) could return incorrectly formatted results when nsim>1.
\item get_intervals(, which="respVar") (and equivalent syntaxes) returned an incorrect result when the family link was not 'identity'.
\item Failure of fitting by sparse-precision method a model specified by a 'corrMatrix' argument with more levels than present in the data.
\item ranef() extractor failed on sparse-precision fits with correlated random effects.
\item prediction variance could be wrong for models with random-coefficient terms and fitted by sparse-precision methods.
\item get_cPredVar() and and confint() no longer stop on fixed-effects models.
\item plot.HLfit() could stop in some cases.
\item eval_replicate() could stop in some cases.
\item Non-default optimizer optim() failed for models with random-coefficient terms. Selecting optim() is not recommended anyway. 
\item Re-emergent glitches in how() when the object had been fitted by a "non-vanilla" call of the fitting function.
}
}
\subsection{DOCUMENTATION}{
\itemize{
\item Substantial sections of the documentation have been rewritten, including new documentation items 'method' and 'random-effects', for easier navigation, updating information, and removing obscurities tracing back to early versions of spaMM. 
\item New documentation 'external-libraries'.
}
}
}
\section{CHANGES UP TO VERSION 3.2.0}{ 
\subsection{NEW FEATURES}{
\itemize{
\item Modified criterion for default selection of sparse precision methods. 
\item is_separated.formula() to test separation with the same arguments as for a fitting function call.
\item Newly exported and documented function dopar() for parallel execution of a function on columns of a matrix. Not logically restricted to mixed-effect applications, hence it can be used more widely.
\item Finer control of verbosity in simulate.HLfit() and predict.HLfit().
\item New 'format' argument for how.HLfit().
}
}
\subsection{NEW DEPENDENCIES}{
\itemize{
\item 'ROI' and 'ROI.plugin.glpk' replace 'lpSolveAPI'.   
}
}
\subsection{FIXED BUGS}{
\itemize{
\item Several errors in simulate() called on COMPoisson-reponse models.
\item Anecdotal glitches in how() when the object was fitted by a "non-vanilla" call of the fitting function.
\item non-default 'REMLformula', as illustrated in example("update.HLfit"), could not be fitted by sparse-precision algorithms.
\item Sparse-precision fits of CAR models could crash for some format of the adjacency matrix.
}
}
}
\section{CHANGES UP TO VERSION 3.1.27}{
\subsection{FIXED BUGS}{
\itemize{
\item simulate(<Gamma-response fit>) simulated according to an incorrect residual variance.
}
}
\subsection{NEW FEATURES}{
\itemize{
\item Faster code for some autoregressive models (mostly for autoregressive LMMs fitted by fitme(), as in tests described at the end of the ''gentle Introduction'' which are ca. 50\% faster).
\item Faster and more memory-efficient code computation of prediction variance for large designs. 'variances$as_tcrossfac_list' control for more memory-efficient (albeit less human-readable) reporting of large prediction covariance matrices.
\item The obsolete functions designL.from.Corr() and sym_eigen() have been removed.
}
}
\subsection{NEW DEPENDENCIES}{
\itemize{
\item 'RSpectra' suggested for a step in fitting autoregressive models.   
}
}
}
\section{CHANGES UP TO VERSION 3.1.2}{
\subsection{NEW FEATURES: Numerical robustness}{
\itemize{
\item Better numerical precision of some computations for binomial models (including a more accurate alternative to binomial(.)$variance(mu)) may avoid numerical imprecisions (or worse) for extreme 'mu'.
\item More rigorous variant of the Levenberg-Marquardt-like algorithm will give better (albeit slower) fits in difficult cases such as binary models with strong spatial correlations, and Gamma GLMMs with response values approaching 1.
\item Default optimizer is not always nloptr(); nloptr() or bobyqa() are now selected depending on input (though nloptr() remains the most common choice).
\item Stricter checking of arguments for adjacency models (the adjacency matrix and any fixed 'rho' parameter).
\item Prediction variances more accurate when correlation matrices of random effecs are nearly singular, for random-coefficient and Matern/Cauchy types of random effects. 
\item Algorithms for fitting COMPoisson-response models made more numerically stable for extreme values of the linear predictor and for low COMP_nu parameter.
}
}
\subsection{NEW FEATURES: Bootstrap computations}{
\itemize{
\item New function eval_replicate() is the new default 'simuland' argument of the likelihood-ratio testing functions, which themselves better handle imperfections of the bootstrap results.
\item Multiple changes to spaMM_boot(), to make it consistently usable in programming, to control parallel::makeCluster(), etc.
\item The '...' argument of LRT() is now passed to spaMM_boot(). 
}
}
\subsection{OTHER NEW FEATURES}{
\itemize{
\item New extractor get_inits_from_fit().
\item By default, spaMM.options() now checks whether in arguments are preexisting options or not. This helps to avoid errors due to typos. The argument 'warn' controls this check.
\item New argument 'verbose' for how().
\item plot_effects() now produces more suitable output than before for binomial models.
\item simulate.HLfit() can perform simulation from linear predictors randomized in various ways (see Details section of its documentation).
\item New argument in get_predVar(., which) provides an alternative extraction method for various components of prediction variance.
\item spaMM.options(separation_max=.) can now be used to control also the checks of separation using the lpSolveAPI package.
}
}
\subsection{NEW FEATURE BY NEW DOCUMENTATION}{
\itemize{
\item 'help("post-fit")' explains how to use multiple-comparison procedures implemented in multcomp::glht on spaMM's fit objects.
}
}
\subsection{USER-LEVEL CHANGES}{
\itemize{
\item fitme() will select iterative methods, rather than "outer" optimization, in some cases where it did not but where it is faster.
\item Computations using the Levenberg-Marquardt-like algorithm tend to be longer than before (there may be exceptions due to other changes).
\item Changes in formal arguments of LRT() (a.k.a. anova.HLfit()) and fixedLRT(). These should be backward compatible, as important removed arguments are now handled through the '...'.
\item SpatialPolygonsDataFrame objects 'worldcountries' and 'oceanmask' removed from the package. See help("seaMask") for availability.
}
}
\subsection{NEW DEPENDENCIES}{
\itemize{
\item 'crayon'; one may not yet see any effect unless verbose=c(TRACE=TRUE) is used.   
\item 'minqa' and 'gmp' moved from Suggests to Imports.   
\item Updated dependency to lpSolveAPI >= 5.5.2.0.17.6 as it seems faster than previous versions on "big problems".
}
}
\subsection{FIXED BUGS}{
\itemize{
\item Interval computation by predict() was not correct for zero-truncated response families.
\item Fitting a truncated negbin model could fail due to numerical precision issues.   
\item HLfit() stopped on non-gaussian response models with several random effects including at least one random-coefficient one (although fitme() allowed to fit such models). Gaussian-response HLfit()'s with the same random-effect structure were imprecise, with typically third-decimal imperfections on the logLik.
}
}
}
\section{CHANGES UP TO VERSION 3.0.0}{
\subsection{NEW FEATURES}{
\itemize{
\item Substantially faster fits of random-coefficient models.
\item Faster procedure for some simple LMMs.
\item New VarCorr() extractor.
}
}
\subsection{FIXED BUGS}{
\itemize{
\item confint() failed (stop) on models fitted by fast procedures suitable for a class of LMMs.
\item The names appearing in certain outputs of lambda values were incorrect.
\item Fits could fail (stop) in the uncommon case where some fixed-effect coefficients were fixed by the user and non-standard REML, as defined by a non-default REML formula, was requested.
\item Prediction variance failed (stop) given more than two random effects.
}
}
\subsection{SUBSTANTIAL CHANGES TO DOCUMENTATION}{
\itemize{
\item 'help(fixed)' completed ('corrPars' information missing), corrected and clarified.
}
}
}
\section{CHANGES UP TO VERSION 2.7.0}{
\subsection{NEW FEATURES}{
\itemize{
\item Trendy stuff: geostatistical analyses based on a Markov Random Field on a lattice are now implemented through the new random-effect syntax 'IMRF' (Interpolated Markov Random Field).
\item Improved implementation of sparse-precision algorithms lead to increases in speed that can exceed 50 \% for fits that use such methods. 
\item New syntax introduced to fit some heteroscedastic non-Gaussian random effects (see help("spaMM")).
\item Better choice of default optimization strategy by fitme() may lead to much faster fits of Gamma-response models.
\item Restriction of [filled.]mapMM() to plots with a single spatial random effect removed, as long as the spatial effects remain together in 2D space.
\item Restriction against predict.HLfit(.,newdata=...), on models with corrMatrix or adjacency terms, lifted when the 'newdata' do not involve levels absent from the correlation structure for such terms.
\item New functions pdep_effects() and plot_effects() to quantify and plot the effects of a variable.  
\item Better control of damping factor of Levenberg-Marquardt algorithm.
\item Better assessment of convergence of Levenberg-Marquardt algorithm for PQL fits.
\item New features of LRT() and spaMM_boot() to facilitate debugging of bootstrap simulations.
\item New package option 'Gamma_min_y' to control numerical problems with low response values in Gamma-response models.
\item update_resp() now handles additional arguments as update() does.
}
}
\subsection{USER-LEVEL CHANGES}{
\itemize{
\item get_boot_response() is deprecated (and removed), and the 'boot_fn' argument of fixedLRT() is deprecated, and its formerly possible value, the function .eval_boot_replicates(), is deprecated (and removed). 
}
}
\subsection{FIXED BUGS}{
\itemize{
\item Predictions function could fail (stop) for some combinations of fixed vs. estimated dispersion parameters.
\item Prediction variance computation could fail (stop) when 'etaFix' argument was used (cryptic case indeed).
\item LRT() could fail (stop) when prior weights were used.
\item anova.HLfit() now properly detect non-nestedness of fixed effects.
\item Information about random number generator was missing from bootstrap output, contrary to the documentation.
\item Incorrect dfs for AIC computation in the case of inner-optimized rho parameter of a CAR model.
}
}
\subsection{SUBSTANTIAL CHANGES TO DOCUMENTATION}{
\itemize{
\item Revised 'help(spaMM_boot)', in particular with respect to the 'simuland' argument.
}
}
}
\section{CHANGES UP TO VERSION 2.6.1}{ 
\subsection{USER-LEVEL CHANGES}{
\itemize{
\item Informative message when plot.HLfit does not produce a 'mean' plot. 
}
}
\subsection{FIXED BUGS}{
\itemize{
\item Fix problem in corMatern example (only on CRAN tests for debian+gcc, which started with a system upgrade to glibc 2.28, and which may reveal weaknesses of nlme and/or nlminb rather than spaMM).
}
}
}
\section{CHANGES UP TO VERSION 2.6.0}{ 
\subsection{NEW FEATURES}{
\itemize{
\item New syntax for autocorrelated random-coefficient effects, allowing e.g. the formula terms Matern(male|.) + Matern(female|.) to fit distinct Matern effects for males and females. See help(spaMM) for details.
\item New function update_resp() to refit on a updated response vector.
\item New simulate() argument 're.form' for finer control of conditioning.
\item New predict() argument 'type' to allow returning the linear predictor.
\item New summary() argument 'verbose' may be convenient when using summary() as an extractor.
\item New function response() to extract the response from the fit.
\item The double-vertical syntax for random effects, ( rhs || lhs ), is now handled. See help(spaMM) for details.
\item HLfit can now fit random-coefficient models by sparse-precision methods.
}
}
\subsection{USER-LEVEL CHANGES}{
\itemize{
\item Faster implementation of some computations in the sparse-precision code.
\item The warnings that occurred when 'e1071' availability was tested, but the package was not available, has been suppressed.
\item Modified approach to controlling numerical precision of prediction covariance matrices, should impact these matrices when they are nearly singular.
\item as_precision() return type documented, and type of the matrix within is itself modified.
}
}
\subsection{FIXED BUGS}{
\itemize{
\item Bug, reintroduced in the previous version, in predict()'ing from models with poly() terms.
\item Bug in simulate(.,type="marginal", newdata=<non-NULL>)
\item confint() had weaknesses, more likely to impact intervals for simple LMMs since the introduction of new fitting algorithms in version 2.5.0.
\item Bug in predict()'ing with an unusual combination of random effects and newdata (documented in test/testthat/test-more-predict.R).
\item Marginal AIC computation did not exclude fixed correlation parameters from number of degrees of freedom.
\item Sparse-precision code could fail when all fixed-effect coefficents were set by the user using the 'etaFix' argument.
\item Failure to catch some incorrect or inefficient calls, notably of HLfit to fit random-coefficient models by sparse-precision methods (which is now a correct call!).
}
}
}
\section{CHANGES UP TO VERSION 2.5.11}{ 
\subsection{NEW FEATURES}{
\itemize{
\item Function get_matrix() to extract various matrices from the fit.
\item Function remove_from_parlist() for manipulating a structured list of parameters.
\item spaMM_boot() now rbind's into a data frame the results of the 'simuland' calls when each call returns a 1-row data frame (the previous behaviour was undocumented and effectively inconsistent accross different parallel backends). 
\item More explicit error message when variables for prior weights are not in the 'data'.
\item logLik(., which="logL_Lap") to compute Laplace approximation using observed Hessian (currently only for family=Gamma(log), for demonstration purposes).
\item Low response variance could result in poor fits. A warning is now issued and several steps have been taken for more robust results.
\item New family() extractor.
}
}
\subsection{FIXED BUGS}{
\itemize{
\item predict() did not generally handle correctly offset terms.
\item A hack perfomed by spaMM on the value of 'options(error=<.>)' could conflict in some cases with the user setting 'options(error=<expression>)'.  
\item Prediction variance computation with a 'fix_X_ZAC.object' failed when the model included a non-autocorrelated random effet.
\item Prior-weights information was not shown in model summaries in some cases where it should.
\item Documentation bug in help("fixed") (with respect to 'etaFix' vs. 'fixed' elements).
\item predict() with 'newdata' did not work from a fit with non-null 'etaFix$beta'.
}
}
}
\section{CHANGES UP TO VERSION 2.5.0}{ 
\subsection{NEW FEATURES}{
\itemize{
\item A faster-fitting procedure is implemented (through fitme()) for fitting a class of LMMs (without heteroscedasticity, without fixed random-effect variances...).
\item The restriction, to models with a single random effect, of the computation of the component of prediction variance for uncertainty in dispersion parameters, has been removed.
\item simulate.HLfit() has new arguments for controlling prior weights and simulation from fits with structured dispersion.
\item EarthChord distance: see help("EarthChord").
\item Models with 'geostatistical' random-effect terms (Matern, Cauchy) can now be fitted by sparse-precision methods (although this is not necessarily very efficient).  
\item Possible improvements in numerical accuracy, notably by sparse-precision methods when (1) fitting models with AR1 terms, that are numerically unstable for extreme AR1 correlation values; (2) in REML fits with low residual variance.
\item Improvement in numerical accuracy of computation of prediction variance when residual variance is low or spatial correlation matrices are nearly-singular.
\item New function MSFDR() implements the multiple-stage false discovery rate method of Benjamini and Gavrilov 2009.
\item Improved .LRT() function (called through fixedLRT()) which makes it more suitable for programming, and which better handles generic-optimization issues.
\item fixedLRT(), LRT() and anova() now call by default spaMM_boot() for the bootstrap. They handle new arguments to control the spaMM_boot() call, allowing finer control of poorly- or non-informative bootstrap samples.
\item is_separated() convenience function to check for (quasi-)separation in binomial-response models (this was previously implemented but not exported nor documented). 
\item Possibly improved initial value for lambda may affect speed and results of some fits.
\item New extractors formula(), terms(), nobs(), and extract.AIC(), being methods for 'HLfit' objects of generics defined by 'stats'. 
\item The syntax '.' can now be used in formulas, e.g. HLfit(y ~ (.)^2, ...)
\item Automatic conversion of correlation matrices to precision factor matrices when sparse-precision algorithm is requested; new as_precision() function for such conversion.
\item More transparent control of the residual variance parameter of a residual dispersion model. However, this affects the default results if the residual dispersion model is a mixed-effect one. Use 'fixed=list(phi=NA)' in the 'resid.model' (forcing estimation of its phi) to recover past results. 
\item More efficient computation of conditional AIC. Users can also control whether it should be computed along with other information criteria. 
\item New extractor how() to cleanly extract information about how a model fit was obtained.
\item An hidden memory-dependent constraint on the size of data that could be fit by fixed-effect models by ML has been alleviated.
\item Numerous improvements for residual dispersion models with random effects makes fitting of the joint residual and main-response models more likely to converge in reasonable time when both models contain spatial terms. 
}
}
\subsection{USER-LEVEL CHANGES}{
\itemize{
\item The Predictor() function is deprecated. See help("Predictor") for details.
\item The designL.from.Corr() function is superseded by mat_sqrt().
\item The structure of the return value of sym_eigen() is modified.
\item 'resid.formula' argument of fitting functions is formally deprecated (it has been declared obsolete since v1.8.0).
\item spaMM now enforces variables for prior weights to be in the 'data'.
}
}
\subsection{PREVIOUSLY DOCUMENTED, BUT NEVERTHELESS MISLEADING FEATURES}{
\itemize{
\item The computation of prediction variance (specifically its component for uncertainty in dispersion parameters) could be incorrect when several random effects were considered. Fixed as noted in 'New Features'.
}
}
\subsection{FIXED BUGS}{
\itemize{
\item The 'random-slope' (Intercept+slope) code was OK, but the more general 'random-coefficients' code (>2 covarying random effects) was not. 
\item Bug in initialization of adjacency correlation parameter, which prevented at least some bootstrap computations. 
\item Fixed a (previously documented) issue in spaMM's update() method, stemming from the distinct interpretation of formulas without explicit fixed-effect terms by spaMM vs. stats. Formula updates should still be carefully checked, as getting them fluent is not on the priority list.
\item Fits failed when spaMM.options(sparse_precision=TRUE) but a corrMatrix was not specified as a precision matrix.
\item Prediction failed from models involving poly(., raw=FALSE) terms.
\item Prediction should use contrast information contained in the fitted object rather than in the session's options().
\item confint() for fixed-effect models has various problems. 
\item An horrendous bug was found in the implementation of Gamma family for inverse link. Fortunately, no reasonable people would use Gamma(inverse)...
}
}
}
\section{CHANGES UP TO VERSION 2.4.34}{
\subsection{NEW FEATURES}{
\itemize{
\item Stricter control of numerical accuracy of COMPoisson-family fits, which does not affect high-nu fits (which were already accurate), but which improve accuracy particularly around nu=1.
\item Use of the 'lpSolveAPI' package is reintroduced to check (quasi-)separation.
\item New extractor get_ranPars() (with limited functionality so far).
\item Parallel computation of some bootstrap computations by foreach+doSNOW has been reinstated. 
\item plot.HLfit has a new argument 'ask', and no longer has a side effect on subsequent plots. 
}
}
\subsection{FIXED BUGS}{
\itemize{
\item Bug in COMPoisson-family fits, which impacted most clearly the results when some response values were zero.
\item Crashing bug could occurs in zero-truncated models (for various cases of the Laplace approximation).
\item Fixed some bootstrap code which could fail (because the original try() statement had disappeared through mis-editing :-( ).
\item simulate(., type = "(ranef|response)") failed, or returned incorrectly formatted results.
}
}
}
\section{CHANGES UP TO VERSION 2.4.8}{
\subsection{NEW FEATURES}{
\itemize{
\item Non-trivial re.form argument now works for models including a 'corrMatrix' term.
\item Precision matrices with more levels than there are in the data are now handled.
}
}
\subsection{FIXED BUGS}{
\itemize{
\item COMPoisson-family fits could crash in some cases (bug introduced in version 2.4.0 in code for initial values).
\item Calls not requesting any form of least-squares fit, such as fitme(y ~ 0, data = <.>), failed. 
\item Prediction interval computation failed for zero-truncated distributions.
\item Prediction variance computation could fail for large newdata ('slicing' algorithm).
}
}
}
\section{CHANGES UP TO VERSION 2.4.0}{
\subsection{NEW FEATURES}{
\itemize{
\item Cauchy spatial correlation family has been implemented.
\item For binomial family the response can now be specified as a factor.
\item Better default selection of methods for autoregressive fits with many random-effect levels may lead to much faster fits.
\item fitme() now selects estimation methods for phi and lambda by default according to nrow(data), which may change details of the results. To overcome the new defaults, fitme() now allows better control of estimation method and of initial values for phi and lambda. 
\item get_predCov_var_fix() may now handle a non-trivial 're.form'.
\item New 'max.print' argument for summary.HLfit().
\item New function spaMM_boot() for parametric bootstrap. 
\item New 'resp_testfn' argument for simulate.HLfit(), to simulate samples that satisfy a given condition.
}
}
\subsection{USER-LEVEL CHANGES}{
\itemize{
\item Revised summary of random-effect parameters.
\item revised format of ranef() screen output for random-coefficient terms, using new method 'print.ranef'.
\item Imports package 'pbapply' to draw progress bars, including in parallel computations; no longer Suggests the package 'foreach' (and optionally 'doSNOW') for the same effect.
}
}
\subsection{FIXED BUGS}{
\itemize{
\item An error was produced when model formulas included a nested random effect (.|./.) followed by a correlated random effect such as Matern(.).  
\item Incorrect prediction variance with non-trivial 're.form' with random-coefficient term(s) or permutations of the random effects.
\item Errors when attempting to fit models with correlated random effects and 'negbin' response family.
\item A wrong algorithm could sometimes be selected for 'inner' estimation of the correlation parameter of CAR models.
\item Wrong documentation for ranef().
\item anova() failed when a fit had non-estimable coefficients because of collinearity.
}
}
\subsection{SUBSTANTIAL CHANGES TO DOCUMENTATION}{
\itemize{
\item Revised 'help(ranFix)' to clarify the usage of argument for fixing parameters in the different fitting functions.
}
}
\subsection{OTHER CHANGES}{
\itemize{
\item Revised default initial values of lambda may slightly change results.
\item Revised default 'fittingFunction' for fixedLRT() may slightly change results, or may call for a user-provided explicit value of 'fittingFunction'.
\item Corrected the poorly written documentation for COMPoisson.
}
}
}
\section{CHANGES UP TO VERSION 2.3.0}{
\subsection{NEW FEATURES}{
\itemize{
\item Zero-truncated variants of Poisson and negative-binomial families are now implemented.
\item Models with first-order autoregressive (AR1) random-effects terms are now fitted.
\item fitme() now fits models with random-coefficient terms, using by default a different strategy than the other fitting functions.
\item Prediction variance computations extended to most models handled by spaMM, including random-slope models.
\item The 'ranFix' argument of HLfit(), and similar 'ranPars' or 'fixed' arguments of other functions, 
      allow a finer control of variance parameters of random effect than before.
\item predict.HLfit() should be substantially faster in most practical cases, particularly when the 'newdata' have many rows.
\item Better approximations are implemented for the infinite sums required for fitting COMPoisson models, in particular for low 'nu'.
\item New extractor get_rankinfo().
\item Some gains in speed, particularly for fit by sparse-precision algorithms, which also benefit from better memory management for large model matrices. 
\item Prediction variance computation now take into account uncertainty in the random effect variance for adjacency random effects. 
}
}
\subsection{USER-LEVEL CHANGES}{
\itemize{
\item Control of optimizers has been revised, with distinct 'optimizer1D' and 'optimizer' in the spaMM.options(); lower case is now used for 'optimizer' keyword in all user-level code; and minqa::bobyqa() can now be called. See ?spaMM.options for further information.
\item Argument 'ZALMatrix' of predictor() is deprecated.
\item Argument 'as_matrix' of get_ZALMatrix() is deprecated.
\item predict.HLfit() gets a new argument 'blockSize'.
\item predict.HLfit( . , binding=NA) newly defined to return a vector instead of a one-column matrix. Attributes will differ too.  
\item Removed dependency on package 'ff'.
}
}
\subsection{FIXED BUGS}{
\itemize{
\item Since version 2.0.0, the ranef() extractor did not return "correlated" effects except for random-coefficient models.
\item Prediction variance computation failed for models involving an adjacency random effect and other random effects. 
}
}
\subsection{SUBSTANTIAL CHANGES TO DOCUMENTATION}{
\itemize{
\item Revision of 'help(autoregressive)' following introduction of AR1 models.
\item Revision of 'help(fitme)' to account for its new features.
\item New 'help(fixed)' gathering information on how to fix parameter values.
\item Simplification and updating of documentation pages for fitting functions.
\item New or updated documentation for 'negbin', 'Tnegbin', 'Poisson', 'Tpoisson' following implementation of the truncated variants.
}
}
\subsection{OTHER CHANGES}{
\itemize{
\item Reimplementation of corMatern (forced by effectively backward-incompatible changes in Rcpp) might very slightly change results.
\item Added '...' argument to LRT(). The function does not yet make anything of it.
}
}
}
\section{CHANGES UP TO VERSION 2.2.0}{
\subsection{NEW FEATURES}{
Widespread improvement in speed for large data sets, thanks to some of the  following changes:
\itemize{
\item Sparse precision matrix algorithm extended to all HGLMs whose linear predictor involves only a restricted set of random effect terms (excluding mainly Matern and random-coefficient terms). It is used by default for conditional autoregressive models, or when the new argument 'covStruct' is used to provide a sparse precision matrix: see help(covStruct). See also help(pedigree) for a type of applications where the latter syntax may be of particular interest, and help(sparse_precision) for more information.
\item Significantly faster algorithm for fitting models with large dense correlation matrices (as in Matern correlation models). 
\item new extractor function get_boot_replicate() to get the response vector for a given bootstrap replicate
\item ‘spaMM_glm()’ and ‘spaMM_glm.fit()’ get the same ‘singular.ok=TRUE’ argument that ‘glm()’ and ‘glm.fit()’ got in R-devel r73380 (2017-09-28).
\item More reliable version of Levenberg-Marquardt algorithm for mixed models (see 'help("LevenbergM")'). Revised rules for using it, so that users typically do not need to think about it.   
}
}
\subsection{USER-LEVEL CHANGES}{
\itemize{
\item fixedLRT() results might be modified in particular when the likelihood has local maxima.
\item Changes in the way the generic optimizers are used may lead to changes in the fitme() and corrHLfit() results.
\item Usage of 'control.corrHLfit' argument has been modified.
\item Fitted values always get their names from the row names of the data.
}
}
\subsection{CORRECTED BUGS}{
\itemize{
\item Fixed numerical issues for binomial(cloglog) family with large values of linear predictor.
\item Random-coefficient models could not be properly fitted by Levenberg-Marquardt algorithm.
}
}
}
\section{CHANGES UP TO VERSION 2.1.6}{
\subsection{NEW FEATURES}{
\itemize{
\item simulate.HLfit() can now simulate samples from a fitted GLMM under the distribution of random effects conditional on the data and given the fixed-effect estimates (using a Gaussian approximation in the GLMM case).
}
}
\subsection{CORRECTED BUGS}{
\itemize{
\item predict.HLfit() could no longer return covariance matrices as a result of a bug introduced in version 2.1.0.
\item 'respVar' computation by predict.HLfit() did not take into account non-identity links, although the documentation suggested it did. 
}
}
}
\section{CHANGES UP TO VERSION 2.1.0}{
\subsection{NEW FEATURES}{
\itemize{
%\item 'HLmethod="SEM"' can now be used to fit probit models if the 'probitgem' package is installed.
\item fitme() can now use a fast sparse matrix algorithm for CAR (adjacency) (G)LMMs. See revised help(adjacency) for further information.
\item confint.HLfit() now computes confidence intervals for GLMs. 
\item spaMM now handles cases where nloptr::nloptr() can be trapped in an infinite loop.
\item residuals.HLfit() has been redefined to bring it closer to stats::residuals.glm(). Up to now it only returned response residuals, now it returns (unscaled) deviance residuals.
\item 'method' argument of spaMM_glm() redefined to allow more experiments (default behaviour is unchanged).
\item simulate.HLfit() now also simulates GLMs.
\item predict(<formula ith random-coefficient term>, <newdata>) is now handled.
\item binomial(cauchit) is now handled.
\item simulate() member function added to COMPoisson(). 
\item simulate() member function modified in negbin() and in spaMM's version of Gamma() to return matrices when nsim>1, as described in help(family).
\item simulate.HLfit() now simulates COMPoisson and negbin (mixed) models.
}
}
\subsection{CORRECTED BUGS}{
\itemize{
\item HLCor() no longer fitted adjacency models in Version 2.0.0.
\item spaMM.options(<new values>) did not return old values.
\item For GLMs, confint.HLfit() (computed no interval but) did return an empty interval.
\item COMPoisson$linkfun sometimes returned theta=0, causing failure of fit.
\item Several problems corrected in degrees of freedom computation for information criteria retruned by get_any_IC(). 
\item Several problems corrected in prediction variance calculation with non-null 're.form'. 
}
}
}
\section{CHANGES UP TO VERSION 2.0.0}{
\subsection{NEW FEATURES}{
\itemize{
\item Core functions have been extensively rewritten to manage expensive matrix computations more efficiently. As a result, version 2.0.0 is 2-3 times faster when correlation matrices of random effects are dense (as typically the case for Matern correlation matrices), and can be much faster when correlation matrices are sparse. 
\item Simple interface for using orthodromic distances: see help(Earth).
\item bootstrap computations can now use parallelisation; see new argument 'nb_cores' for fixedLRT() and related functions. 
\item The term for uncertainty in dispersion parameters is now included in computation of prediction variance (1) by default from objects returned by fitme() (special default 'control' argument no longer needed); and (2) for models with multiple random effects.
\item New extractors get_predCov_var_fix() and preprocess_fix_corr() for more efficient computation of specific blocks of prediction covariance matrix. 
\item get_ZALMatrix() gets a new argument 'as_matrix' to control the return value.
\item confint() now works for HLCor fits.
\item nloptr() call replaces optim() call in random-slope code, with notable improvements in some toy examples. Random-slope models with >1 predictor variables, such as (x+y|z), are now fitted.
\item New function stripHLfit() to remove some large objects from a fit object before saving it to file.
\item Time to fit a model is now recorded in the fit object (as element 'fit_time').
}
}
\subsection{USER-LEVEL CHANGES}{
\itemize{
\item Structure of bootstrap output in return value of fixedLRT() modified. 
\item Structure of ranef() result modified. 
\item Several efforts to reduce the size of fitted objects.
\item Likelihood-ratio chi-square statistic renamed 'chi2_LR' (was 'LR2') in return value of anova() and some bootstrap outputs. 
\item Raster.Rd removed to reduced package dependencies. A link to a web page is instead provided. 
\item spaMM now provides a slightly less informative diagnostic for (quasi-)separation in binary response models, but no longer 'imports' package lpSolveAPI, instead 'suggesting' e1071.   
\item Tracing stuff, such as 'trace' argument of some functions, has been removed.
\item Possible values of 'control.corrHLfit$Optimizer' are (as before) 'L-BFGS-B' and (new) 'NLOPT_LN_BOBYQA'. nlminb can no longer be selected.  
\item get_ZALMatrix() is now a standard function rather than part of the object resulting from a fit.
\item fitme() now fits CAR model parameters using its own optimization procedure rather than through HLfit's procedures. 
\item Message added to warn users about the interpretation of formulas without fixed effects.
\item More meaningful display of results for random-coefficients models; new argument 'details' for summary.HLfit().
\item ML fit summaries no longer display a 'p_bv' value, to prevent confusion about what this value means. 
}
}
\subsection{BUGS}{
\itemize{
\item Version 1.10.0 failed on R < 3.2.0 because cBind() and rBind() calls from earlier spaMM versions have been replaced by cbind() and rbind() calls (see help(cBind)). The Depends field is now corrected to R (>= 3.2.0).  
\item A few stopping bugs in predict() have been corrected.  
\item LRT() did not work on HLCor fit objects. This is corrected.
}
}
}
\section{CHANGES UP TO VERSION 1.10.0}{
\subsection{NEW FEATURES}{
\itemize{
\item spaMM now fits models with random effects in the residual dispersion model, and thus implements a class of \dQuote{double hierarchical GLMs}.
\item fitme() now fits models with offset in the dispersion model.
\item Modified convergence criteria for the iterative algorithm in HLfit().
\item Examples run generally faster by about 30\%. 
\item New extractor function get_intervals() and get_ZALMatrix().
\item Intervals returned by predict.HLfit() (and hence by get_intervals()) use Student's t distribution when the fitted model is a simple linear model.
\item Prediction variance can now be computed from fits of models with a corrMatern term.
\item An HLfit object's $APHLs$p_bv for a *fixed-effect* model now returns the restricted likelihood (for estimation of residual dispersion parameters), which it has never done so far (out of disregard for this case). This allows REML fits of fixed-effect models using fitme() (but has no bearing on similar fits using HLfit()). 
\item More informative summary() when non-default spatial distance function, or non-default prior weights, are used.
\item New 'refit' control in fitme().
\item New extractor function get_RLRTSim_args() for use with the RLRsim package to perform small-sample LRT of the presence of a random effect.
}
}
\subsection{USER-LEVEL CHANGES}{
\itemize{
\item Summary() outlook revised, in particular to more accurately display the meaning of phi for Gamma response. 
\item The objects returned by the fitting functions should be of much smaller size than before. 
\item Better reporting of terms that have to be ignored in the computation of prediction variance.  
\item New 'inverse.Gamma', 'spaMM-conventions' and 'member-extractor' documentation pages. HLfit() documentation revised.  
}
}
\subsection{BUGS}{
\itemize{
\item Some numerical issues that could affect GLM fitting with HLfit() have been fixed.
\item More accurate detection of unidentifiable random effect components, when residual variance formulas involve a pure offset term.
\item confint() on a LMM fitted by fitme() was returning incorrect results. This is corrected.
}
}
}
\section{CHANGES UP TO VERSION 1.9.16}{
\subsection{NEW FEATURES}{
\itemize{
\item Computation of prediction variance has been reimplemented so that it requires storing smaller matrices than before. Further, some terms required for this computation are computed only when prediction variance is requested, rather than being always returned by HLfit(). This allows some gains in computation time for large data sets. 
\item simulate.HLfit() gets new arguments and prints some information, in anticipation of future developments.
\item spaMM.colors() gets new argument 'adjustcolor_args' for finer control.
\item Filled plot functions gets new argument 'plot.title'.
}
}
\subsection{BUGS}{
\itemize{
\item spaMM_glm(), which has been designed to 'rescue' stats::glm() when the latter fails, did not rescue it for some Gamma(log) GLMs. This is corrected.
\item Computation of AIC (and other information criteria) failed in many cases. This is corrected.
\item confint() did not correctly handle models with several random effects. This is corrected.
\item fitme() does not fit random-slope models, nor models with offset in the dispersion model, but this was not documented. This is corrected.
\item C++ code has been fixed for compatibility with recent standards.
}
}
}
%
\section{CHANGES UP TO VERSION 1.9.0}{
\subsection{NEW FEATURES}{
\itemize{
\item New fitme() function, intended to become the main fitting function in spaMM, can be substantially faster than HLfit() or corrHLfit(), particularly for large data sets, but can occasionally be slower.
\item Some other moderate but noticeable improvements in speed. 
\item Users no longer need to explicitly request AIC computation when calling the fit of a model. The new extractor function get_any_IC() (or its more conventionnally named alias, AIC()) can be used to compute them from a fitted model object.
\item Negative binomial response models can be fitted by 'family=negbin()'.
\item The syntax 'family=COMPoisson()', with unknown 'nu' parameter, works. This means that COMPoisson() models can now be fitted with unknown values of this parameter, by corrHLfit() or by the new fitme() function. The same is true for the parameter of the 'negbin' family. 
\item Initial or fixed values of the variance of random effects can be controlled independently when there are several random effects.
\item spaMM_glm() gains a new argument 'strict' for finer control. spaMM_glm.fit() has increased default 'control$maxit'.  
\item filled.mapMM() now returns a predictor of the response surface. It also gains a new argument 'smoothObject' that takes such a predictor, allowing for faster redrawing.   
\item fixedLRT() gets a new argument 'fittingFunction' to select between corrHLfit() and fitme().
}
}
\subsection{USER-LEVEL CHANGES}{
\itemize{
\item Revised usage of 'corrMatrix' argument of HLCor().
}
}
\subsection{BUGS}{
\itemize{
\item 'spaMM', as well as 'RcppEigen', had become broken for Windows 32bit. See https://github.com/RcppCore/RcppEigen/issues/34 for details. A suitable compilation flag is now used to overcome this issue.   
\item simulate.HLfit() could fail when a design matrix was stored as a _M_atrix. This is corrected.
}
}
}
\section{CHANGES UP TO VERSION 1.8.0}{
\subsection{NEW FEATURES}{
\itemize{
% negative.binomial
%\item Stochastic Expectation-Maximization (SEM) algorithm for binary probit model (to be documented)
\item predict.HLfit() can return prediction intervals.
\item New spaMM_glm() and spaMM_glm.fit() functions provide a robust alternative (using the Levenberg-Marquardt algorithm) to glm.fit() when the latter diverges. They also circumvent the initial-value issues of glm().
\item Slight general increase in speed. 
\item The Conway-Maxwell-Poisson family is implemented.
\item New extractor functions get_fixefVar(), get_predVar(), get_residVar() and get_respVar().
\item Subsetting the result of a predict() call also subsets its variance attributes, thanks to a new internal 'Extract' method `[.predictions`.
\item spaMM is protected from disturbing numerical problems of the Gamma()$dev.resids() function, thanks to an internal redefinition of the Gamma family.
\item New 'control.glm' argument of HLfit(), to control calls to spaMM_glm() within HLfit(). 
}
}
\subsection{USER-LEVEL CHANGES}{
\itemize{
\item Once-again revamped 'variances' argument and documentation of predict.HLfit(). The 'variances' argument of filled.mapMM() is replaced by 'variance'.
\item New argument 'resid.model' replaces 'resid.formula' and 'resid.family' specifications in the main fitting function.
\item The default behaviour of mapMM() for argument 'map.asp=NULL' has been modified.
\item The obscure control parameters 'RHOMAX' and 'NUMAX' previously accessible through spaMM.options() are no longer needed and have been removed.
\item Better control of some plot margins.
\item Fitting random-slope models may be substantially less slow (but still slow).
}
}
\subsection{BUGS}{
\itemize{
\item Printing bug for non-standard REML formula in summary.HLfit() corrected. 
\item Conditions where predict() return an object of class 'predictions' corrected to avoid a printing bug. 
\item The 'xrange' and 'yrange' arguments of spaMMplot2D() were not effective. This is corrected.
\item Printing bug in summary.fixedLRT() corrected.
\item spaMM.options(maxLambda=. . .) did not fully control extreme lambda estimates. The revised code is more effective, and also better bypasses some numerical problems that can occur for such values.
\item Fitting phi using the corrHLfit() outer optimization mechanism was inefficiently implemented and could return some confusing output. This is corrected.  
}
}
}
\section{CHANGES UP TO VERSION 1.7.2}{
\subsection{NEW FEATURES}{
\itemize{
\item Less inefficient code for fixed-effect models, faster and handling larger data sets.
\item Changed default value of 'upper' argument of optim call() within corrHLfit(), and better reporting of issues with this call.
}
}
\subsection{USER-LEVEL CHANGES}{
\itemize{
\item Argument 'which' replaces arguments 'REML' in logLik().
\item New control 'break_conv_logL' in argument control.HLfit of function HLfit(). Mainly for development purposes.
\item Improved reporting of convergence issues in HLfit().
\item More comfortable print method for vcov() value.
}
}
\subsection{BUGS}{
\itemize{
\item This new version includes several improvements of predict(), the most notable being a correction of prediction variances for new values of predictor variables, which were incorrect in CRAN version 1.6.2 (but not in CRAN version 1.5.1). Many other more trivial bugs or features of predict() have been corrected.   
\item summary.fixedLRT() failed to display correctly the bootstrap results. This is corrected.
\item A bug in C++ code, introduced in version  1.7.1, has been corrected in version 1.7.2; and code applying some Eigen functions on a matrix with a dimension equal to zero (generating an UBSAN runtime error but a correct retrun value in other conditions) has been fixed.
}
}
}
\section{CHANGES UP TO VERSION 1.6.2}{
\subsection{NEW FEATURES}{
\itemize{
\item Prediction variance can now account for uncertainty is dispersion parameters (though not for all models; see Details of 'disp' in the 'predict' documentation).
\item Now allows prediction for a large number of spatial positions by storing some large correlation matrices on disk rather than in memory (using the ff package, and as controlled by the spaMM option 'ff_threshold'). (NB: this was removed in a later version)
\item Convenience extractor functions 'dev_resids' and 'deviance' have been added. 
\item New data objects 'worldcountries' and 'oceanmask' added for better plots of spatial predictions (see the new documentation '?raster' for plot of spatial predictions using rasters and geographical projections). (PS: removed since, see version 3.1.2) 
\item More efficient implementation of a conditional autoregressive (CAR) model. % (including in SEM algorithm) (to be documented)
\item 're.form' argument of 'predict' function has been extended to handle values of class 'formula'.
\item etaFix$beta usage has been extensively revised to provide a convenient interface for fixing (some of) the fixed effects (beta coefficients) values.
\item If anyone on Earth ever tried inverse.Gamma("-1/mu") random effects, s-he might have been confronted with negative estimates of u_h. This problem is inherent to the algorithms described in the HGLM literature, but spaMM now brings a solution, which allows one to fit more securely models with v ~ Gamma( lambda/(lambda+1), lambda+1). A new link identity is also allowed for Gamma random effects [i.e., v=u ~ Gamma(1/lambda,lambda) ], albeit with the restriction that lambda is < 1.
\item New 'symSVD' argument for designL.from.Corr function.
\item New keyword 'predict' for 'which' argument of 'plot.HLfit', for easier plotting of prediction diagnostic plot.
\item 'make.scaled.dist' has been renamed 'make_scaled_dist', and gets a new 'return_matrix' argument for programming purposes.
}
}
\subsection{USER-LEVEL CHANGES}{
\itemize{
\item Substantial gains in speed for moderately large data sets (hundreds of lines), in particular by minimizing the use of the QR module from Eigen.
\item Users no longer have to think about pre-computing prediction coefficients. 'predictionCoeffs' function has been removed.
      'coeffs' argument of predict.HLfit() has been removed. 
\item More compact printing of predict() results.
\item The output format of information criterias computed by HLfit, and related documentation, have been revised.
}
}
\subsection{BUGS}{
\itemize{
\item confint() did not correctly handle PQL fits. This is corrected. In addition, it better handles some non-convergence issues.  
\item Gamma response models with non-unit 'prior.weights' produced an inaccurate correction from EQL to 'HL(1,1)'. This is corrected.
\item anova() and LRT() declared some nested models as non-nested, and sometimes stopped. This is corrected.
\item HLfit could stop on some fits involving moderately large sparse design matrices (Z) for random effects. This is corrected.
\item A printing bug for ranef(), which prevented the display of its result when lambda was fixed, has been corrected 
}
}
}
%
\section{CHANGES IN VERSION 1.5.1}{
\subsection{USER-LEVEL CHANGES}{
\itemize{
\item New example_maxtime package option for controlling execution of examples.
}
}
\subsection{TESTS}{
\itemize{
\item Tests are now controlled by a local environment variable _LOCAL_TESTS_ so that they are not run on CRAN.
}
}
}
%
\section{CHANGES IN VERSION 1.5.0}{
\subsection{NEW FEATURES}{\itemize{
\item Faster code for (some) sparse models with a single random effect on large data sets. It is also possible (by using \code{spaMM.options(QRmethod=Matrix::qr)}) to use sparse matrix algorithms in other models, but this will be slow in many applications.
Various other changes were made to speed up computation in dense models (by about 10\% in the tests).

\item The non-Euclidean distance functions implemented in the \code{proxy} package (e.g. Geodesic) can now be used (although only Geodesic was tested, and the Matern correlation function with arbitrary rho and nu is not mathematically feasible for all distance functions). See \code{control.dist} argument of \code{HLCor} and \code{corrHLfit}.
}}
\subsection{USER-LEVEL CHANGES}{\itemize{
\item Reorganized graphic functions, with new \code{spaMMplot2D} function and new \code{envir} argument for \code{mapMM.HLfit}.
\item \code{predict.HLfit} modified for better consistency with other widely used predict methods: added \code{newdata} argument and modified default behaviour of \code{binding} argument. 
\item New \code{getDistMat} extractor for extracting the distance matrix used in the Matérn model 
\item This NEWS.Rd file is created
}}
\subsection{BUG FIXES}{\itemize{
\item Correction of a bug that caused incorrect prediction on new data for 1D (only) spatial model. 
}}
}