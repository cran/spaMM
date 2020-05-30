\name{external-libraries}
\alias{external-libraries}
\title{
Installing external libraries
}
\description{
spaMM is conceived to minimize installation issues but it nevertheless suggests using some external libraries. These are all accessed through R packages so their installation should be easy when installing binary packages.  The Details below give hints for installing packages from source. For all cases not considered below, help yourself.
}
\details{
The \code{ROI.plugin.glpk} package requires the \code{Rglpk} package, which itself requires the external \code{glpk} library. For the latter, Debian-ists and alikes should \code{sudo apt-get install libglpk-dev}. MacOSX users should \code{brew install glpk} if using \code{brew}; Windows users may need to install glpk from \url{https://sourceforge.net/projects/winglpk/}.

The \code{nloptr} package requires the external \code{NLopt} library. If you need to install it from source on Windows and depending on the version of g++ included in your installed version of the Rtools, you may need to be a geek, as available precompiled static libraries may not link. 
%\url{https://github.com/aadler/nloptr/blob/Windows_2.6_plus/INSTALL_Windows.md}.
}