\name{spaMM-internal}

% TODO % * retry outer beta estim ? Without gradient... hum. see 'outerbetaML' test
% * Something about .predictions methods when predict returns a data frame ?

\alias{gridIMRF} % __F I X M E___
%\alias{.setCluster} % __F I X M E___ document ? It is exported for Infusion

\alias{projpath} % for local tests

\alias{glm.nodev.fit} % _F I X M E__ special version for COMPoisson

\alias{fitme_body}
\alias{fitmv_body}
\alias{HLfit.obj}
\alias{.HLfit_body_augZXy} % alias bc exported bc it must be exported to be traceable
\alias{.solve_IRLS_as_ZX} %  alias bc exported bc it must be exported to be traceable
\alias{.solve_IRLS_as_spprec} % alias bc exported bc it must be exported to be traceable
\alias{.is_spprec} % exported bc widely used in tests

% More aliases of exported internal functions because R CMD check (debian at least) otherwise interprets them as undocumented user-level:

\alias{.NB_shapeFn}
\alias{.NB_shapeInv}
\alias{.ZWZt}
\alias{.calc_CARdispGammaGLM}
\alias{.calc_dispGammaGLM}
\alias{.crossprod}
\alias{.dispFn}
\alias{.dispInv}
\alias{.evalWrapper}
\alias{.is_spprec_fit}
\alias{.modify_list}
\alias{.nuFn}
\alias{.nuInv}
\alias{.rhoFn}
\alias{.rhoInv}
\alias{.safe_opt}
\alias{.setCluster}
\alias{.stripOffset}
\alias{.stripRanefs}
\alias{.tcrossprod}
\alias{.unloads4spaMM}
\alias{.wrap_Ltri_t_chol}

\alias{dimnames.bigq} % hidden method for base:: generic

% exported from nlme
\alias{Initialize} 
\alias{Variogram} 
\alias{coef} 
\alias{coef<-} 
\alias{corFactor} 
\alias{corMatrix} 
\alias{getCovariate} 
\alias{logDet} 
\alias{recalc} 
\alias{spaMM_Gamma} 
\alias{get_HLCorcall} % fixme: deserves to be documented as of independent interest (cf exampl of usage in probitgem)
\alias{HLCor.obj} 
%\alias{ULI} % cf def for what to do
\alias{corrHLfit_body}
\alias{corrMM.LRT} % in Ecography code, redefined as a stop()
\alias{getCall.HLfit} 
\alias{model.offset.HLfit} % not an S3 method
\alias{getPar}
\alias{HLCor_body}
\alias{HLfit_body}
\alias{HLfit_body_old}
\alias{makeTicks} % used by probitgem and Infusion but not spaMM itself. Either I define it in each of the two other packages or I keep it here
\alias{Matern.corr} % back compat, should be deleted later
\alias{niceLabels} % used by blackbox
\alias{overcat} % used in probitgem. Infusion and blackbox both have .overcat()

\alias{crack} % not doc'ed but used in an example :-)
\alias{recond} % not doc'ed and used in test

\title{Internal spaMM Functions}
\description{
  Internal spaMM functions
}
\details{
  These are not to be called by the user, or are waiting for documentation to be written.
}
\keyword{ internal }
