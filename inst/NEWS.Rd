\name{NEWS}
\title{spaMM News}
\encoding{UTF-8}
% cf ?news for the format and usage of this file
%
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
\item New extractor function get_RLRTSim_args() for use with the LRTSim package to perform small-sample LRT of the presence of a random effect.
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
\item The syntax 'family=COMPoisson()', with unknown 'nu' parameter, works. This means that COMpoisson() models can now be fitted with unknown values of this parameter, by corrHLfit() or by the new fitme() function. The same is true for the parameter of the 'negbin' family. 
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
\item Now allows prediction for a large number of spatial positions by storing some large correlation matrices on disk rather than in memory (using the ff package, and as controlled by the spaMM option 'ff_threshold').
\item Convenience extractor functions 'dev_resids' and 'deviance' have been added. 
\item New data objects 'worldcountries' and 'oceanmask' added for better plots of spatial predictions (see the new documentation '?raster' for plot of spatial predictions using rasters and geographical projections). 
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