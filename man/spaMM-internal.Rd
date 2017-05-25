\name{spaMM-internal}

\alias{generateName} 
\alias{makenewname}

\alias{merge_parlist}
\alias{modify_list}

\alias{fitme_body}
\alias{calc_maxrange}
\alias{HLfit.obj}




%\alias{\%id*id\%}
%\alias{as_call_family}
\alias{Bartlett.robust} % not currently used 
\alias{Beta} 
\alias{calc.w.resid}
\alias{calc_wAugX}
\alias{calc.plot.dims}
\alias{calcD2hDv2}
\alias{calcPHI}
\alias{calc_dispGammaGLM}
\alias{calcRanefPars}
\alias{calcUniqueGeo}
\alias{checkRandLinkS}
\alias{checkRandLink}
%\alias{checkRespFam}
%\alias{compare.model.structures}
\alias{corr.notEQL.lambda}
%\alias{Dim.corMatern} % not currently used
%\alias{eval.update.call} % not currently used
\alias{eval.corrEst.args}
\alias{provide_init_lambdas}
\alias{spaMM_Gamma} 
%\alias{generateFileName} % only a stop() for obsolescence
\alias{getHLfit}
\alias{get_valid_beta_coefs}
\alias{get_HLCorcall}
\alias{HLCor.obj} 
\alias{HL_process.args} 
\alias{HL_strictify} 
\alias{HLframes}
\alias{Initialize.corMatern} % probably augments nlme 
\alias{initialize_v_h}
\alias{ULI}
\alias{Variogram.corMatern} % probably augments nlme 
\alias{asStandardFormula} 
\alias{calcNewCorrs}
\alias{calcPredVar}
\alias{calcResidVar}
\alias{calc_logdisp_cov}
\alias{canonizeRanPars} 
\alias{clean_cs} 
\alias{checkDistMatrix} 
\alias{coef.corMatern} % probably augments nlme
\alias{coef<-.corMatern} % probably augments nlme
\alias{corFactor.corMatern} % probably augments nlme
\alias{corMatrix.corMatern} % probably augments nlme
\alias{corrHLfit_body}
\alias{corrMM.LRT} 
\alias{d2mudeta2fn} 
\alias{dispFn}
%\alias{calc_etaGLMblob}
\alias{calc_CARdispGammaGLM} 
\alias{dispInv} 
%\alias{eval_gain_LevM_GLM}
\alias{extract.check.coords} 
\alias{findOffset} 
\alias{findSpatial}
\alias{parseBars} 
\alias{findbarsMM} 
\alias{findplus} % not currently used (except recursively by itself)
\alias{fullrho}
\alias{generateInitPhi}
\alias{getCall.HLfit} 
\alias{get_beta_cov_any_version}
\alias{getCovariate.corMatern} % probably augments nlme
\alias{getgroups} % not currently used
\alias{getPar}
\alias{getProcessed}
\alias{getValidData}
%\alias{get_methods_ranefs}
\alias{get_ZAlist}
\alias{get_LSmatrix}
\alias{glm.reformat}
\alias{HLCor_body}
\alias{HLfit_body}
\alias{is.identity} 
\alias{largeLambdaMessages} 
\alias{legend_lambda}
\alias{LevenbergMstepCallingCpp}
\alias{locoptim}
\alias{logDet.corMatern} % probably augments nlme
\alias{LogAbsDetWrap}
\alias{makeCheckGeoMatrices}
\alias{makeCovEst1} 
\alias{makeCovEst2} 
\alias{makeLowerUpper}
\alias{makeTicks}
\alias{MaternCorr.ff}
\alias{Matern.corr} % back compat, should be deleted later
\alias{Matrix_times_Dvec}
\alias{muetafn}
\alias{nextPoints}
\alias{niceLabels}
\alias{nobarsMM} 
\alias{nobarsMM_} 
\alias{nobarsNooffset} 
\alias{noNonSpatialbarsMM}
\alias{noOffset} 
\alias{noOffset_} 
\alias{.noRanef} 
\alias{nuFn} 
\alias{nuInv} 
\alias{overcat}
\alias{pastefrom} 
\alias{plotScale}
\alias{post_process_family}
\alias{preprocess}
\alias{print.arglist}
\alias{print.predictor}
\alias{print.vcov.HLfit}
\alias{print.predictions}
\alias{process_resglm_list}
\alias{recalc.corMatern} 
%\alias{reformat_verbose}
\alias{resize_lambda}
\alias{rhoFn} 
\alias{rhoInv}
\alias{safesolve.qr.matrix} 
\alias{safesolve.qr.vector} 
\alias{selectLoglfn}
%\alias{pMVN}
%\alias{separator} 
\alias{setProcessed}
\alias{setControl}
\alias{Sigwrapper}
\alias{Sig_times_b}
\alias{singularSigmaMessagesStop}
\alias{slashTerms} 
\alias{solveWrap.matrix}
\alias{solveWrap.vector}
\alias{internal_LRT} % not renamed .LRT as foreach does not handle it
\alias{subbarsMM} 
%\alias{theta.mu.canonical}
\alias{theta.mu.conjugate} 
\alias{tryCatch.W.E}
\alias{validData}

\alias{crack} % not doc'ed but used in an example :-)


\title{Internal spaMM Functions}
\description{
  Internal spaMM functions
}
\details{
  These are not to be called by the user, or are waiting for documentation to be written.
}
\keyword{ internal }
