// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
//
// Zigg.h:  Base class for different Ziggurart implementations
//
// This file is taken from RcppZiggurat ( Copyright (C) 2013  Dirk Eddelbuettel).
//
// RcppZiggurat is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// RcppZiggurat is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with RcppZiggurat.  If not, see <http://www.gnu.org/licenses/>.

#ifndef RcppZiggurat__Zigg_h
#define RcppZiggurat__Zigg_h

#include <cmath>
#include <stdlib.h>  // FR: RcppZiggurat does not need this... !
#include <stdint.h>             // not allowed to use cstdint as it needs C++11

namespace Ziggurat {

    class Zigg {
    public:
        virtual ~Zigg() {};
        virtual void setSeed(const uint32_t s) = 0;
        // no getSeed() as GSL has none
        virtual double norm() = 0;
    };
}
    
#endif
