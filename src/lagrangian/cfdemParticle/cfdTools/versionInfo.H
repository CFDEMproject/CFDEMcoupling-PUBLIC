/*---------------------------------------------------------------------------*\
    CFDEMcoupling - Open Source CFD-DEM coupling

    CFDEMcoupling is part of the CFDEMproject
    www.cfdem.com
                                Christoph Goniva, christoph.goniva@cfdem.com
                                Copyright 2009-2012 JKU Linz
                                Copyright 2012-     DCS Computing GmbH, Linz
-------------------------------------------------------------------------------
License
    This file is part of CFDEMcoupling.

    CFDEMcoupling is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    CFDEMcoupling is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with CFDEMcoupling; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Description
    This code is designed to realize coupled CFD-DEM simulations using LIGGGHTS
    and OpenFOAM(R). Note: this code is not part of OpenFOAM(R) (see DISCLAIMER).

\*---------------------------------------------------------------------------*/

#ifndef versionInfo_H
#define versionInfo_H

word CFDEMversion="cfdem-3.8.1";
word compatibleLIGGGHTSversion="3.8.0";
word OFversion="5.x-commit-538044ac05c4672b37c7df607dca1116fa88df88";

Info << "\nCFDEMcoupling version: " << CFDEMversion <<  endl;
Info << ", compatible to LIGGGHTS version: " << compatibleLIGGGHTSversion <<  endl;
Info << ", compatible to OF version and build: " << OFversion <<  endl;

#endif
