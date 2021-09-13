\name{spaMM-internal}

% TODO ___F I X M E___
% * implement outer optime of beta starting from the etaFix (+ REMLformula) code => avoid LevM costs for COMP ? 
% * Something about .predictions methods when predict returns a data frame ?

\alias{projpath} % for local tests

\alias{glm.nodev.fit} % __F I X M E__ possible temporary

\alias{fitme_body}
\alias{fitmv_body}
\alias{HLfit.obj}
\alias{.HLfit_body_augZXy} % alias bc exported bc it must be exported to be traceable
\alias{.solve_IRLS_as_ZX} %  alias bc exported bc it must be exported to be traceable
\alias{.solve_IRLS_as_spprec} % alias bc exported bc it must be exported to be traceable

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
\alias{getPar}
\alias{HLCor_body}
\alias{HLfit_body}
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
