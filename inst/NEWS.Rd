\name{NEWS}
\title{spaMM News}
\encoding{UTF-8}
% cf ?news for the format and usage of this file
%
\section{CHANGES UP TO VERSION 4.4.16}{  
\subsection{NEW FEATURES}{
\itemize{
\item New predict.HLfit(., na.action) argument to control handling of missing predictors in the 'newdata'.
\item Some computations in predict.HLfit() may be performed in parallel, as controlled by new argument 'cluster_args'.
\item mapMM() and spaMM.filled.contour() now return some graphic elements of the plot.
\item Code made compatible with anticipated change in result of NCOL(NULL) in base R. 
}
}
\subsection{FIXED BUGS}{
\itemize{
\item Several issues with post-fit functions on LMs with sparse design matrix (as may happen by default when the predictor variable is a factor with many levels).
\item Internal default value of 'map.asp' in spaMMplot2D() was not always appropriate, resulting in a failure of the plot.
\item pdep_effects() and plot_effects() did not handle factors well enough, and did not handle fits with non-estimable fixed-effect coefficients.
\item Some fits could fail with recent versions of the Matrix package.
\item Some prediction varaince computation could fail for AR1 models.
}
}
}
\section{CHANGES UP TO VERSION 4.4.0}{  
\subsection{NEW FEATURES}{
\itemize{
\item predict(<multivariate-response fit>, re.form=<formula>, ...) and simulate(<idem>) now more thoroughly implemented. 
\item Easier control of multivariate-response bootstraps made possible by several changes. Notably, simulate4boot() gets a '...' argument, for additional arguments now passed to simulate.HLfit(). See also help("respName").
\item Corr() return value is modified for some models, notably IMRF ones. 
\item Improved choice of sparse versus dense matrix methods for models including several Matern( 0+<indicator for level of factor> | <coords>) terms (and surely also for other correlation models than Matern).
\item Normalization of 'IMRF' models (including 'MaternIMRFa') on irregular grids, so that the random effect is homoscedastic, is now possible.
\item New DoF() extractor for degrees of freedom.
}
}
\subsection{FIXED BUGS}{
\itemize{
\item Several bugs have been fixed for simulate(., newdata=.) for multivariate-response fits. 
\item Fixed a stopping bug due to a typo, in rarely executed code handling a case where the Cholesky factorisation of a model matrix failed.
\item Lower alpha value for 'MaternIMRFa' terms set to 1+1e-08, as inla::inla.spde2.matern() is not numerically stable in alpha=1. 
\item Multivariate-response fits mixing random effects implemented through corrFamily objects together with random effects not implemented in this way could fail (stop). 
\item A workaround has been implemented to handle a bug in geometry::tsearch (possibly related to https://github.com/davidcsterratt/geometry/issues/46 and /57) which had the effect that some data located within a mesh defining an IMRF term were treated as if out of the mesh.
}
}
\subsection{DEPENDENCIES}{
\itemize{
\item spaMM now imports the 'backports' package.
}
}
}
\section{CHANGES UP TO VERSION 4.3.20}{ %     
\subsection{NEW FEATURES}{
\itemize{
\item LRT() (a.k.a anova.HLfit() with two input fits) revised to better suggest variants of LR tests suitable for the models compared, and to more accurately (and sometimes rather verbosely) report possible problems.
\item New simulate4boot() function, useful for precomputing bootstrap samples to be used by spaMM_boot(). The latter function gets a new argument 'boot_samples' to input such precomputed samples.
\item New spaMM2boot() function, differing from spaMM_boot() by the format of its return value, which is suitable for input for boot::boot.ci().
\item New simulate_ranef() function.
\item New lev2bool() convenience function for conversion of factor to 0/1 variable.
\item div_info() revised to provide more clear information; but this feature is still experimental.
}
}
\subsection{ANTICIPATING VERSION 1.6.0 OF 'Matrix' PACKAGE}{
\itemize{
\item Backward-incompatible changes will be introduced in Matrix v1.6.0. spaMM code has been preemptively modified to prevent  bugs that would result from these changes, including a silent one (sigh). Still, users might occasionally need to recompute some fit objects created with older versions of Matrix.
}
}
\subsection{BUG FIXES}{
\itemize{
\item Problems with 'composite' Matern random effects: error when a level is found only once in the data, and predict() failures.
\item LRT(<fit of GLM>, test="F"), and probably a few other cases of the ANOVA tables for fixed-effect models, could fail when the design matrix is singular.
}
}
}
\section{CHANGES UP TO VERSION 4.3.0}{ %   
\subsection{NEW FEATURES}{
\itemize{
\item Beta-binomial models can be fitted ('betabin' family).
\item New 'clinics' data set.
\item The residVar() extractor now handles missing information in its 'newdata' argument 
}
}
\subsection{USER-LEVEL CHANGES}{
\itemize{
\item deviance() now returns weighted residuals when relevant, for consistency with stats::deviance. Its documentation now details this distinction, with several examples.
}
}
\subsection{BUG FIXES}{
\itemize{
\item Result of dev_resids() and deviance() extractors was not appropriately documented for non-GLM families. 
\item Residual-variance extractors had not fully taken account of new features of spaMM (such as new response families) for a long time and could fail for some combinations of arguments (e.g., get_residVar() with 'newdata' for a multivariate-response fit with non-GLM families). 
\item summary.HLfit() failed on fixed-effect models with non-GLM responses families and fixed (by 'etaFix' argument) fixed-effect coefficients.
\item A fit could fail (stop) when the model formula included a factor with many levels as fixed effect, and a random effect with dense correlation matrix. 
}
}
\subsection{DEPENDENCIES}{
\itemize{
\item spaMM now depends on R>=3.5.0.
}
}
}
\section{CHANGES UP TO VERSION 4.2.0}{  
\subsection{NEW FEATURES}{
\itemize{
\item New fixed-effect terms defined by functions X.GCA() and X.antisym() for dyadic interactions. 
\item Residual dispersion models can be fitted for the negbin1, negbin2 and beta response family. 
\item Most restrictions on composite random-effect terms have been removed. For example, composite AR1 or ARp random effects can now be fitted. 
\item fitmv() gets new argument 'X2X' allowing to fit shared fixed-effect coefficients among different submodels of a multivariate-response fit. 
\item LRT() now handles comparisons of models with any form of residual-dispersion model. 
\item (Hopefully) improved control of starting values for fixed-effect coefficients, relevant for some difficult cases with log link.   
\item numInfo() more correctly handles partially-fixed random-coefficient parameters.
\item anova() when used as interface to lmerTest methods, more correctly handles partially-fixed random-coefficient parameters.
\item numInfo() gets new 'sing' argument controlling the computation and display of small eigenvalues of the information matrix.
\item Extended handling, by predict() and simulate() for multivariate-response models, of missing information in 'newdata'.
\item New get_fittedPars() extractor.
\item The range of response values is now checked for beta response family.
\item New combinepar() interface for parallel computations, always calling foreach to handle any 'control$.combine' argument.
\item Faster COMPoisson-response fits (from a few to more the 50 percent depending on arguments).
\item New 'Leuca' data set.
}
}
\subsection{USER-LEVEL CHANGES}{
\itemize{
\item Some names in output of prediction intervals and of confidence intervals depend on confidence levels, whose conversion to string depends on the 'OutDec' global option. The relevant spaMM code now locally ensures that OutDec is '.' irrespective of any globally-set value, to favour portability of spaMM results.    
\item The package option 'CMP_asympto_cond' should now be a function rather than a (quoted) expression. Formally, this is a backward-incompatible change as the options is documented (though changing it was not encouraged).   
}
}
\subsection{FIXED BUGS}{
\itemize{
\item Binomial GLMs with non-canonical link could no longer be fitted (! -- mixed-effect models were OK). 
\item For some response family, fixed-effect models could not, or no longer, handle design matrices in sparse Matrix format. 
\item A computation was incorrect for COMPoisson distribution with huge expectation. Fixing it has had no visible effect so far. 
\item Fits of simple mixed-effect models with random-coefficient terms were subject to two potential problems:\cr
\verb{ * } When fitted by "inner" estimation methods [in particular, by HLfit() or HLCor() rather than fitme()], the fit could stop in cases where the observed information matrix differs from the expected one (e.g, for the Gamma(log) response family). \cr
\verb{ * } Using fitme() to fit such models worked for all response families, but a profile confint() procedure for fixed-effect coefficients on the resulting fit object could then stop.\cr \cr
Although using HLfit() to fit such models now works for all response families, the confint() procedure does not work on the resulting fit object. This problem is cleanly reported, with a suggestion to fit the model by fitme(). 
\item All problems with deviance computations for the 'negbin1' family had not been solved in the previous version, more specifically for the zero-truncated family, where it has been corrected in this version (numerical inaccuracies can still occur in the case of jointly vanishing 'scale' and 'mu'). Even more than previously, the imperfections of these computations may be unnoticeable in practice. 
}
}
\subsection{DEPENDENCIES}{
\itemize{
\item 'memoise' suggested dependency has been removed.
}
}
}
\section{CHANGES UP TO VERSION 4.1.20}{ 
\subsection{NEW FEATURES}{
\itemize{
\item LRT() can now perform comparisons of fits differing only by their residual-dispersion models, when the latter are fixed-effect models.
\item Change in the default return value of the fixef() extractor for rank-deficient mixed-effect models, creating a potential incompatibility with some code using this extractor, but actually helpful for the next feature:
\item Various improvements of the procedure for anova() tables, in particular (1) handling of models with rank-deficient model matrices; (2) handling of multivariate-response models; (3) handling of various possible problems in computation or use of numerical information matrix (LMM case).
\item numInfo() is more robust, in particular by checking whether some parameter estimates are at boundary of their possible range (see Details of the function documentation). 
}
}
\subsection{USER-LEVEL CHANGES}{
\itemize{
\item 'check' argument of drop1.HLfit() renamed as 'check_marg'.
}
}
\subsection{FIXED BUGS}{
\itemize{
\item A change in Matrix v1.5-3 has exposed a bug in the code for sparse-precision methods, which would result in incorrect or even failed (stop) fits. spaMM results with previous versions of Matrix were correct.
\item A 'new feature' in Matrix v1.5-2 ("Products of diagonalMatrix with [RT]sparseMatrix now preserve storage...") also broke a procedure for IMRF models.
\item Requesting prediction variances with the (awkward, but formally correct) argument 'variances=list(linPred=FALSE, cov=TRUE)' generated nonsense.  
\item predict(., <\emph{large newdata (> 2000 rows)}>, binding= < \emph{some name} >, ...) retained only the first 2000 element of any attribute for variances or intervals. 
}
}
}
\section{CHANGES UP TO VERSION 4.1.0}{ 
\subsection{NEW FEATURES}{
\itemize{
\item drop1.HLfit method implemented, to return tables of comparisons of a model with single-term deletions from it, with specific procedures for LMs, GLMs and LMMs. The LMM case is handled by interfacing the Satterthwaite method in \code{lmerTest}.
\item New function numInfo() to compute the information matrix for (ideally) all parameters for mean response model, using numerical derivatives. Experimental.
\item New antisym(<.>|<.>) random effect to fit anti-symmetric dyadic interactions, as in Bradley-Terry models.
\item Pseudo-R2 computation has been implemented, see help("pseudoR2"). Please don't (ab)use it.
}
}
\subsection{USER-LEVEL CHANGES}{
\itemize{
\item 'fixed' argument now available, to fix model parameters, in all fitting functions; older 'ranFix' and 'ranPars' arguments of some functions are still available.
}
}
\subsection{FIXED BUGS}{
\itemize{
\item Deviance computation for the 'negbin1' family was broken, though the effect could be unnoticeable in fits of mixed-effect models (where this was only used to provide some initial values, and not even for all 'negbin1' fits). 
}
}
}
\section{CHANGES UP TO VERSION 4.0.0}{ 
\subsection{NEW FEATURES}{
\itemize{
\item Laplace approximation using observed, rather than expected, information matrix (see help("method")) is no longer "experimental". It is now available for all response families, including in multivariate-response fits, and is now the default method, changing the results of fits using response family with non-canonical link (hence the new "major" version number).
\item anova.HLfit() can now return traditional ANOVA tables, with specific procedures for LMs, GLMs and LMMs. The LMM case is handled by interfacing the Satterthwaite method in \code{lmerTest::anova}. The whole functionality is experimental and unstable as later versions might considered alternative default procedures.
\item Access to procedures implemented in the lmerTest package (contrast tests with effective degrees of freedom, and anova tables) is also possible through new function as_LMLT() allowing conversion of a spaMM fit object to input suitable for these procedures. 
\item Prediction for levels of a 'corrMatrix' random effect that were not present in the fitted data is now possible (if these levels are present in the input corrMatrix). Previously this was not possible because only the subset of the corrMatrix for levels present in the fitted data was kept in the fitted model object.   
\item Additional extractors have been implemented ('df.residual' method) or extended ('weights', 'residuals').
}
}
\subsection{USER-LEVEL CHANGES}{
\itemize{
\item 'negbin2' is the new preferred names of the negative-binomial GLM family (to better distinguish it from the non-GLM 'negbin1' family). 'negbin' and 'Tnegbin' are now aliases for negbin2.  
\item logLik(., which="logL_Lap") is deprecated.
\item revised get_ranPars(), including new which="fitted" option.
}
}
\subsection{BUG FIXES}{
\itemize{
\item A bug was introduced in version 3.12.0 affecting the fits of nested random effects using the (. | . \%in\% .) syntax. The examples in help("GxE") no longer returned correct results.
\item Fit results may be modified for some non-GLM response family models by sparse-correlation methods (although so far the fixed bug affected only the progress of the fit but not its final result, except in a case where the information matrix is not positive definite at the attained maximum, so that the fit result is suspect anyway).
\item the get_ranPars() extractor did not handle well random-coefficient parameters
\item Fitting beta family models could fail when some response value was exactly 1/2.  
}
}
}
\section{CHANGES UP TO VERSION 3.13.0}{ 
\subsection{NEW FEATURES}{
\itemize{
\item Non-GLM response families can now be fitted, the first two implemented being the 'negbin1' and 'beta_resp'. This is experimental and are limitations (notably, this is not yet possible in multivariate-response fits).
\item plot_effects() now handles the case where the focal variable is a factor.
\item AIC.HLfit() now handles several fitted-model objects as described for the stats::AIC generic (the AIC.HLfit formals have been modified for that purpose but no backward-incompatibility should result). 
\item New 'which' value "fixef_left_ginv" handled by get_matrix() extractor.
}
}
\subsection{DOCUMENTATION}{
\itemize{
\item Added example using ggplot2 in help("plot_effects").
\item The documentation previously failed to mention summary.HLfit() as a useful *extractor*.
}
}
\subsection{USER-LEVEL CHANGES}{
\itemize{
\item For truncated GLM families (Tpoisson, Tnegbin), the definition of deviance residuals has been changed [for details see help("Tpoisson")]. This changes the deviances reported by a glm() fit using these families; and it may theoretically affect mixed-effect model fits, as the deviance may affect initial values of some estimates, although no effect on fit results has been observed.
\item Prediction intervals (as opposed to the response intervals) will be on the linear predictor or response scale depending on \code{type}. Previously they were always on response scale.
\item get_ranPars(., which=NULL) (the default) did not return all lambda estimates (as get_ranPars(., which="lambda") does), but only "outer-estimated ones" (as get_ranPars(., which="outer_lambda") does). This poorly documented feature has been modified.   
}
}
\subsection{FIXED BUGS}{
\itemize{
\item Degrees of freedom (as used by LRT(), notably) could be incorrectly counted for models mixing "outer-estimated" random-coefficient models and other, "inner-estimated", random-effect parameters.
\item The elements of the "table" attribute of the result of confint.HLfit() could be misnomed.
\item predict.HLfit(., re.form = <not NULL>) failed when some fixed-effect coefficients estimates were NA, or were set by etaFix argument.  
\item Non-default specification of the type(s) of confidence interval in a confint.HLfit() call was ignored (not passed to boot::boot.ci()).
}
}
}
\section{CHANGES UP TO VERSION 3.12.0}{ 
\subsection{NEW FEATURES}{
\itemize{
\item Experimental extension of geostatistical random effects to include "composite" Matern terms, meaning in particular Matern(<LHS>|<RHS>) terms where the LHS is a numeric variable or factor.
\item Even more experimental implementation of generic handling of <corrFamily>(.|<RHS>) for nested RHS (see last example in help("corrFamily-design")).
\item Partial implementation of Laplace approximation using observed, rather than expected, information matrix (see help("method"), vastly rewritten on that occasion).
\item hatvalues() extractor can now return results for multivariate-response fits.
\item More robust convergence of algorithm for GLMs fitted by fitme().
}
}
\subsection{FIXED BUGS}{
\itemize{
\item Predictions with 'newdata' for "adjacency" models were mostly incorrect.
\item Nested geostatistical random effects [such as Matern( 1 | . \%in\% . )] were likely to fail except in "balanced" nested designs.
}
}
\subsection{USER-LEVEL CHANGES}{
\itemize{
\item Minor adjustment to pdep_effects(), so that it now returns a single row of output when given a single-row 'newdata' or 'newdata' with a single value of the focal (numeric) variable.
\item Better control of extreme random-effect variances ('lambda') may lead to numerical changes in some fits, mainly those with low variances.  
\item Improved convergence warnings.
}
}
}
\section{CHANGES UP TO VERSION 3.11.14}{ 
\subsection{FIXED COMPATIBILITY ISSUE}{
\itemize{
\item Use of ...names() had created an overlooked dependency on R >= 4.1. 
}
}
\subsection{USER-LEVEL CHANGES}{
\itemize{
\item Much faster fits by COMPoisson() (>10 times generally, and even >25 times faster in some fits with default "loglambda" link).   
}
}
\subsection{FIXED BUGS}{
\itemize{
\item get_inits_from_fit() did not properly work on all residual dispersion models. 
}
}
}
\section{CHANGES UP TO VERSION 3.11.3}{       
\subsection{NEW FEATURES}{
\itemize{
\item The features of the experimental 'corrFamily' random effect have been extended for easier definition of new random-effect correlation models, and this has been used to add predefined random-effect models, as next listed.
\item One can fit autoregressive AR(p) and ARMA(p,q) models of any order, using formula terms of the form ARp(<.>|<.>) or ARMA(<.>|<.>).
\item One can fit IMRF models, including (in contrast to IMRF random effect terms) estimation of the 'alpha' parameter (analogous to the Matern smoothness) of Lindren et al's (2011) approximation of the Matern model by a Gaussian Markov random field, by using the formula term MaternIMRFa(<.>|<.>).   
\item One can fit some models for diallel experiments and other symmetric dyadic interactions, using formula terms of the form diallel(<.>|<.>) or ranGCA(<.>|<.>).
\item The fitting functions and fixedLRT() get a new argument 'weights.form' to specify prior weights in a safe way in all programming contexts. 
}
}
\subsection{USER-LEVEL CHANGES}{
\itemize{
\item Backward-compatible change in the 'verbose' argument of HLfit() and HLCor(): 'verbose=c(inner=<.>)' should now be preferred to 'verbose=c(trace=<.>)', to avoid confusion with the more useful 'verbose=c(TRACE=<.>)' of other fitting functions.   
\item LRT.HLfit() (a.k.a. anova.HLfit()) will return full results in more cases where it is unable to check that models are nested (with proper messages).
\item The computation of a large matrix is now avoided in the summary of sparse-preision fits. 
}
}
\subsection{DEPENDENCIES}{
\itemize{
\item New dependency to geometry package [which removes a "soft" dependency to the not-on-CRAN 'INLA' package (function inla.spde.make.A) when fitting IMRFs].
}
}
\subsection{FIXED BUG}{
\itemize{
\item 'prior.weights' in multivariate-response models could not be specified in the way recommended in help("good-practice").    
\item Invalid code for default y-axis label for plot_effects() with binomial model where the LHS of the formula is a single variable rather than a cbind(.,.) expression.
\item The replacement function for stats::gaussian( < non-identity link > )$initialize() introduced in the previous version did not handle the case where fixed-effect coefficients were fully fixed.
}
}
}
\section{CHANGES UP TO VERSION 3.10.0}{       
\subsection{NEW FEATURES}{
\itemize{
\item New experimental 'corrFamily' random effect to fit correlated random effects belonging to a user-defined correlation model.
\item New convenience functions extreme_eig() and regularize() to handle nearly-singular matrices (prompted by 'corrFamily' development).
\item The general-purpose dopar() and dofuture() functions can now execute code depending on other packages than spaMM whatever the parallel backends these functions call.
\item spaMM fitting functions, as well as spaMM_glm.fit(), rested on stats::gaussian()$initialize() to find good starting values for gaussian-response models with log or inverse link, which means that they would sometimes fail, forcing users to provide initial values. This is no longer necessary (see exponential regression example in help("spaMM_glm")).   
}
}
\subsection{USER-LEVEL CHANGES}{
\itemize{
\item Faster fitme/fitmv fits of random-coefficient models constrained to diagonal covariance matrix by a fixed=list(ranCoefs=<....>) argument: useful notably for (<factor>|<RHS>) or (<0+factor>|<RHS>) terms with several levels of the factor.
\item Cryptic changes (to some loop termination condition and to a matrix operation) to improve the fitting time of various slow fits (with a fourfold reduction in one extreme case).
\item Revised procedures for regularization of nearly-singular matrices, affecting mat_sqrt() results for them. This may in principle impact some fits.
}
}
\subsection{BUG FIXES}{
\itemize{
\item confint() could fail in some programming contexts.
\item Fitting some LMMs without fixed effects (not really L'M'Ms then...) could fail.
\item summary() could fail for GLMs without estimated fixed-effect coefficient.
\item The Levenberg-Marquardt-like algorithm could seem stuck in rares cases with large random effects (thousands of levels).
\item The shape of the negbin() family was ignored when the shape argument was a call to some function (and likewise for COMPoisson's nu).
}
}
\subsection{DOCUMENTATION}{
\itemize{
\item Example of analysing interactions by random effects, in new documentation topic help("GxE").
\item Use of the (experimental) dummy() factor specifier from the lme4 package is now documented.
}
}
}
\section{CHANGES UP TO VERSION 3.9.25}{       
\subsection{NEW FEATURES}{
\itemize{
\item New argument filled.mapMM(., return.) to control the return value of this function.
\item New argument 'focal_values' for functions pdef_effects() and plot_effects(), to control values at which predictions are computed.
\item These News for version 3.9.0 said "Much more memory-efficient (and faster) fits for 'nested-geostatistical' random effects..." but one more improvement was needed to bring down memory requirements over all steps of a fit.
}
}
\subsection{USER-LEVEL CHANGES}{
\itemize{
\item The fixed-effect coefficients returned by get_inits_from_fit(from, template) no longer include values for coefficients that were NA in fixef(template).  
}
}
\subsection{BUG FIXES}{
\itemize{
\item Bootstrap procedures in e.g. LRT() could fail when the fixed-effect coefficients of the "full" model fit, as returned by fixef(), included NA's (as expected when some model coefficients cannot be estimated from the data). This has been fixed by the above-mentioned change in the return value of get_inits_from_fit().
}
}
}
\section{CHANGES UP TO VERSION 3.9.13}{     
\subsection{NEW FEATURES}{
\itemize{
\item Faster fits by sparse-precision methods when (1) several random effects are present (2) random-coefficient terms are present. (Very) small speed gains in many other fits.
}
}
\subsection{BUG FIXES (ALMOST)}{
\itemize{
\item Fixed an inconsequential clang compiler warning (... indeed, no real bug fix in this version).  
}
}
\subsection{DEPENDENCIES}{
\itemize{
\item New 'Suggests': 'memoise', memoisation being optionally used on a few spaMM functions.
}
}
}
\section{CHANGES UP TO VERSION 3.9.0}{    
\subsection{NEW FEATURES}{
\itemize{
\item New type of random effect handled: "composite" random effects such as corrMatrix(<variable>|ID) or corrMatrix(mv(1,2)|ID), combining features of a random-coefficient model and of a corrMatrix (see help("composite-ranef")). Such random effects are considered in quantitative genetic analyses of multivariate response, for example.
\item residuals.HLfit() can now also return standardized deviance residuals used in fitting methods and in diagnostic plots.
\item Much more memory-efficient (and faster) fits for 'nested-geostatistical' random effects (e.g., Matern(1|long+lat \%in\% year) ).
\item The 'control.HLfit' argument of fitting functions now handles an 'algebra' parameter for easier control of matrix-algebraic methods. See help("algebra") for possible usage. 
\item Procedures which rest on the stats::gaussian(<"log"|"inverse">) family object to fit models with such response family to data with zero(es) in the response may fail to find valid starting values, which is not totally defensible. spaMM now circumvents this obstruction. Whether asking for fits of such data by such families is a good idea is another matter; the principle adopted here is to allow the fit and let its results contribute to answer this question.
\item pdep_effects() now handles 'length.out=0L' with the effect that predictions are made for actual values of the predictor variables, then allowing for adding informative plots of fit-residuals on the partial-dependence plots. 
}
}
\subsection{USER-LEVEL CHANGES}{
\itemize{
\item For random-coefficient terms fitted by sparse-precision methods, there was a discrepancy between some elements stored in the fitted object. Neither the fits nor predictions were affected, but trying to use together the results of extractors get_ZALMatrix() and ranef.HLfit() would have been confusing. This discrepancy has been removed. 
\item plot.HLfit() should now return a diagnostic plot for standardized deviance residuals in all cases. 
\item More consistent default value of the residual variance of the residual variance model (sic): it is now set to a fixed value in mixed-effect residual-variance models specified by resid.model=<formula>, as was already effectively so for other specifications of the residual-dispersion model. 
}
}
\subsection{FIXED BUGS}{
\itemize{
\item predict(., <with newdata>) results were incorrect for models with multivariate random-coefficient terms (i.e., using the mv() syntax) fitted by sparse-precision methods.
\item Prediction variance was incorrect for models with random-coefficient terms.
\item Prediction variance computation could fail when all random-effect variances were fixed.
\item residVar(., which="var") had incorrect results for gaussian-response models. This did not affect response-variance computations nor get_residVar(), which did not call this extractor. residVar() could also fail on objects with an offset-only formula for phi.
}
}
}
\section{CHANGES UP TO VERSION 3.8.23}{    
\subsection{NEW FEATURES}{
\itemize{
\item Improved condition for selection of sparse-precision methods for models with 'corrMatrix' terms, with quite substantial gains in speed for affected fits. 
\item Numerical approximations for COMPoisson-family models have been substantially revised, with disappointing gains in speed. Additional gains are possible if the 'memoise' package has been installed on your R system. 
\item The 'control.HLfit' argument of fitting functions now handles a 'sparse_precision' parameter (previously controllable only globally through spaMM.options()). See help("sparse_precision") for possible usage. 
\item The criterion for switching to Levenberg-Marquardt-like algorithm has been slightly modified and hopefully improved. 
\item New 'Gryphon' data set.
}
}
\subsection{USER-LEVEL CHANGES}{
\itemize{
\item Long-deprecated argument 'uniqueGeo' of HLCor() is now formally removed.
}
}
\subsection{FIXED BUGS}{
\itemize{
\item Prediction variance computations could fail for models with random-coefficient terms and not fitted by sparse-precision methods,
due to a bug introduced in v3.7.48 (i.e., CRAN versions 3.8.0 and 3.8.9)
\item Better control of arguments to avoid bugs [eval() looking for the IMRF 'model' in the global environment] when fitting is called in a programming context.  
\item Bug (stop) when calling predict.HLfit() with both 'intervals' calculation and binding of point predictions to data.
}
}
\subsection{DOCUMENTATION}{
\itemize{
\item Substantially revised documentation for "sparse_precision".
}
}
}
\section{CHANGES UP TO VERSION 3.8.9}{   
\subsection{NEW FEATURES}{
\itemize{
\item Modified selection between ''outer'' optimization versus ''inner'' iterations to acheive faster default fitme() fits of some non-spatial GLMMs. 
\item The default cluster type is now "PSOCK" in all cases.
\item New interface dofuture() for parallelisation. See help("wrap_parallel") for differences from dopar().
\item The 'control.HLfit' argument of fitting functions now handles a 'sparse_precision' parameter (previously controllable only globally through spaMM.options()). See help("sparse_precision") for possible usage. 
\item Improved condition for selection of sparse-precision methods for models with 'corrMatrix' terms. 
}
}
\subsection{FIXED BUGS}{
\itemize{
\item Typo in a recently introduced check of geographical coordinates led to a warning when abs(longitude) values >90 occurred.  
\item Tests using the 'adjlg' data set no longer worked. It was the code generating the data set, rather than any spaMM procedure, that required fixing.
\item Using fork clusters provided statistically correct results but could not be controlled as documented (the random number generator and progress bar could not be properly controlled, and incorrect warnings could be displayed). 
}
}
\subsection{DEPENDENCIES}{
\itemize{
\item New 'Suggests': 'future' and 'future.apply'.
}
}
}
\section{CHANGES UP TO VERSION 3.8.0}{    
\subsection{NEW FEATURES}{
\itemize{
\item Parallel computations can now use forking via parallel::mclapply.
\item The conditional AIC can now be evaluated by a bootstrap procedure.
\item Checks that coordinates are consistent with longitude and latitude for Earth and EarthChord distances.
\item Previous version had become slower for some spatial GLMMs without residual dispersion parameter. This is fixed.  
}
}
\subsection{FIXED BUGS}{
\itemize{
\item Fits with "adjacency" random effects failed when data were missing for some locations. 
\item predict(., binding= < char string >, blockSize = < any value > ) no longer worked when predicting more than 'blockSize' values. 
\item For Gamma-response model, get_residVar() (and same computation through call to predict()) returned the dispersion parameter phi rather than the residual variances. 
\item Use of the dot in a formula, to refer "to the remaining variables contained in data" (cf ?terms.formula), is now handled, with a warning about its limitations. It could previously lead to errors (stop), e.g. if all such variables were then removed  ( ~ . - ... - ...) to leave only the Intercept.
}
}
\subsection{DOCUMENTATION}{
\itemize{
\item Substantially revised documentation for prediction variances, with new ?predVar page, and for AIC computations.
}
}
}
\section{CHANGES UP TO VERSION 3.7.34}{   
\subsection{NEW FEATURES}{
\itemize{
\item Non-canonical links implemented for the COMPoisson family. 
\item More precise computation of prediction variance for some cases with vanishing residual variance (and possibly a few other cases). 
\item New 'which' value in residVar(.,which="fam_parm") to extract the family dispersion parameter of COMPoisson and [T]negbin families. 
\item A "table" attribute in the format of the stats::confint() result is added to the confint.HLfit() result for profile confidence intervals.
}
}
\subsection{USER-LEVEL CHANGES}{
\itemize{
\item Possibly large speed gains for a small class of REML fits on large data sets/random effects with many levels, due to these being now more likely to be fitted by ``outer'' optimization. 
\item Fits of models with a resid.model of the form ~ 0 + offset(.) should also be faster. This affects the 'mean_fit' performed by IsoriX::isofit(), for example. 
\item how.HLfit()'s message about fitting method is less cryptic. Its return value is unchanged.
\item Warning added about potential back-compatibility issue in some multinomial fits due to changes in default optimization strategy.
}
}
\subsection{FIXED BUGS}{
\itemize{
\item Issue with missing data in multivariate-response fits, that could show up when summary()zing the fit.
\item Prediction variance calculations with a re.form of the form '[no LHS] ~ <single random-effect>' ignored that random-effect (this did not occur with a dummy LHS nor with several random effects).
\item VarCorr() would report NA instead of the variance for random-slope terms without correlation parameter [cf terms of the form ( 0 + X |  <grouping factor> )].
\item summary.HLfit() could fail in some programming contexts.
}
}
}
\section{CHANGES UP TO VERSION 3.7.0}{  
\subsection{NEW FEATURES}{
\itemize{
\item Expressions such as I(y>0) or even y>0 are now accepted in the left-hand side of a model formula (up to now update_resp() and bootstrap functions failed on them, and I(.) was explicitly disallowed). 
}
}
\subsection{USER-LEVEL CHANGES}{
\itemize{
\item 5-10\% speed gains in many fits (at least non-spatial ones), and even more, notably (1) for Levenberg-Marquard-like algorithm (hence for fits affected by convergence issues, which tend to be the longest ones), and (2) for fits by sparse-precision methods for (2a) Gamma(log) GLMMs (>fifthfold in a large problem!) or (2b) other GLMMs when fitted by PQL rather than ML.
\item In addition, the criterion for using ``inner'' iterative methods rather than ``outer'' optimization ones for fitting dispersion parameters has been modified, resulting in substantially faster fits in some non-spatial models. 
\item New diagnostic message pointing a potential issue when fitting models with group-specific geostatistical effects without matching Intercept terms.
}
}
\subsection{FIXED BUGS}{
\itemize{
\item Back-compatibility fixes for summaries and AIC computations on old fit objects (tested against v2.5.0 objects).
}
}
}
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
\item It was given a new scrutiny, leading to a display bug fix, a revised fitted() output, improved summary() information, a new how.HLfitlist() extractor, and revised documentation. But fitmv() might ultimately be more useful to fit multinomial data. 
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
\item Sparse-precision code could fail when all fixed-effect coefficients were set by the user using the 'etaFix' argument.
\item Failure to catch some incorrect or inefficient calls, notably of HLfit to fit random-coefficient models by sparse-precision methods (which is now a correct call!).
}
}
}
\section{CHANGES UP TO VERSION 2.5.11}{ 
\subsection{NEW FEATURES}{
\itemize{
\item Function get_matrix() to extract various matrices from the fit.
\item Function remove_from_parlist() for manipulating a structured list of parameters.
\item spaMM_boot() now rbind's into a data frame the results of the 'simuland' calls when each call returns a 1-row data frame (the previous behaviour was undocumented and effectively inconsistent across different parallel backends). 
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
\item Numerous improvements for residual dispersion models with random effects makes fitting of the joint residual and mean-response models more likely to converge in reasonable time when both models contain spatial terms. 
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