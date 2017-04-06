\name{spaMM.colors}
\alias{spaMM.colors}

\title{
  A flashy color palette.
}
\description{
  \code{spaMM.colors} is the default color palette for some color plots in \code{spaMM}.
}
\usage{
spaMM.colors(n = 64, redshift = 1, adjustcolor_args=NULL)
}
\arguments{
  \item{n}{Number of color levels returned by the function. A calling graphic function with argument \code{nlevels} will typically 
  take the first (i.e., bluest) \code{nlevels} color levels. If \code{n}<\code{nlevels}, the color levels are recycled
  }
  \item{redshift}{The higher it is, the more the palette blushes....}
  \item{adjustcolor_args}{Either NULL or a list of arguments for \code{\link[grDevices]{adjustcolor}}, in which case \code{adjustcolor} is called to modify \code{spaMM.colors}'s default vector of colors. See the documentation of the latter function for further information. All arguments except \code{col} are possible.}
}
\value{
A vector giving the colors in a hexadecimal format.
}
\details{
If you don't like this color palette, have a look at the various ones provided by the \code{fields} package.
}

\examples{
## see mapMM examples
}
%- \keyword{ print }
