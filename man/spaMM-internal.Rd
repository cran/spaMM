\name{spaMM-internal}

%\alias{generateName} 
\alias{makenewname} % fixme can be removed when next (>08/2017) version of Infusion is on CRAN

%\alias{modify_list}

\alias{fitme_body}
\alias{HLfit.obj}


\alias{noNonSpatialbarsMM} % fixme can be removed when next (>08/2017) version of Infusion is on CRAN
\alias{LogAbsDetWrap} % fixme can be removed when next (>08/2017) version of Infusion is on CRAN


%\alias{Bartlett.robust} % not currently used 
\alias{Beta} 
%\alias{calc_wAugX}
%\alias{calcPHI}


%\alias{calc_dispGammaGLM} % used by probitgem !
%\alias{calc_CARdispGammaGLM} % used by probitgem !





%\alias{calcRanefPars}
%\alias{calcUniqueGeo}
%\alias{corr.notEQL.lambda}
%\alias{eval.corrEst.args}
\alias{spaMM_Gamma} 
\alias{get_HLCorcall} % fixme: deserves to be documented as of independent interest (cf exampl of usage in probitgem)
\alias{HLCor.obj} 
%\alias{HL_process.args} 
%\alias{HL_strictify} 
\alias{ULI} % cf def for what to do
%\alias{asStandardFormula} 
%\alias{calcNewCorrs}
%\alias{calcPredVar}
%\alias{calcResidVar}
%\alias{calc_logdisp_cov}
\alias{corrHLfit_body}
\alias{corrMM.LRT} % in Ecography code, redefined as a stop()
%\alias{dispFn}
%\alias{dispInv} 
%\alias{generateInitPhi}
\alias{getCall.HLfit} 
%\alias{get_beta_cov_any_version}
\alias{getPar}
%\alias{getProcessed}
%\alias{get_LSmatrix}
%\alias{glm.reformat}
\alias{HLCor_body}
\alias{HLfit_body}
%\alias{locoptim}
%\alias{makeCovEst1} 
%\alias{makeCovEst2} 
%\alias{makeLowerUpper}
\alias{makeTicks} % used by probitgem and Infusion but not spaMM itself. Either I define it in each of the two other packages or I keep it here
\alias{Matern.corr} % back compat, should be deleted later
\alias{niceLabels} % used by blackbox
%\alias{nobarsNooffset} 
%\alias{noOffset} 
%\alias{noOffset_} 
%\alias{.noRanef} 
%\alias{nuFn} 
%\alias{nuInv} 
\alias{overcat} % in probitgem. Was in Infusion and blackbox (not in the future > 26/08/2017 versions) so cannot be removed yet.
\alias{pastefrom} 
%\alias{plotScale}
%\alias{post_process_family}
%\alias{preprocess}
%\alias{process_resglm_list}
%\alias{rhoFn} 
%\alias{rhoInv}
%\alias{setProcessed}
%\alias{internal_LRT} % ---> .LRT

\alias{crack} % not doc'ed but used in an example :-)


\title{Internal spaMM Functions}
\description{
  Internal spaMM functions
}
\details{
  These are not to be called by the user, or are waiting for documentation to be written.
}
\keyword{ internal }
