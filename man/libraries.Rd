\name{external-libraries}
\alias{external-libraries}
\title{
Installing external libraries
}
\description{
spaMM is conceived to minimize installation issues but it nevertheless suggests using some external libraries. These are all accessed through R packages so their installation should be easy when installing binary packages.  The Details below give hints for installing packages from source. They may all be obsolete if you are using the Rtools42 on Windows. For all cases not considered below, help yourself. If you are using the Rtools40 on Windows, you should have a look at the package manager in the Rtools40 bash shell.
}
\details{
The \code{ROI.plugin.glpk} package requires the \code{Rglpk} package, which itself requires the external \code{glpk} library. For the latter, Debian-ists and alikes should \code{sudo apt-get install libglpk-dev}. MacOSX users should \code{brew install glpk} if using \code{brew};  Windows users should try using\cr
 \code{pacman -S mingw-w64-x86_64-glpk} in the Rtools40 bash shell, together with\cr
 \code{Sys.setenv(GLPK_HOME = "$(MINGW_PREFIX)")} in the R session (but I have not fully tested this; previously I had to install glpk from \url{https://sourceforge.net/projects/winglpk/}).

The \code{nloptr} package requires the external \code{NLopt} library. Windows users should try using\cr
 \code{pacman -S mingw-w64-x86_64-nlopt} in the Rtools40 bash shell (but again I have not fully tested this; see also the README of \code{nloptr}). 
%\url{https://github.com/aadler/nloptr/blob/Windows_2.6_plus/INSTALL_Windows.md}.
}
