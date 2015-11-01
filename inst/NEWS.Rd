\name{NEWS}
\title{spaMM News}
\encoding{UTF-8}
% cf ?news for the format and usage of this file
%
\section{CHANGES UP TO VERSION 1.6}{
\subsection{NEW FEATURES}{
\itemize{
%\item Stochastic Expectation-Maximization (SEM) algorithm for binary probit model (to be documented)
\item Prediction variance can now account for uncertainty is dispersion parameters (though not for all models; see Details of 'dispVar' in the 'predict' documentation).
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