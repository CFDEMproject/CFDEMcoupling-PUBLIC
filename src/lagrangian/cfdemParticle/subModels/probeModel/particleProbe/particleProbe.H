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

Class
    particleProbe

SourceFiles
    particleProbe.C

\*---------------------------------------------------------------------------*/

#ifndef particleProbe_H
#define particleProbe_H

#include "probeModel.H"
#include "fvCFD.H"
#include "polyMesh.H"
#include "cfdemCloud.H"
#include "OFstream.H"
#include <vector>
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class particleProbe Declaration
\*---------------------------------------------------------------------------*/

class particleProbe
:
    public probeModel
{

private:

    dictionary propsDict_;

    mutable word name_;

    cfdemCloud& particleCloud_;

    bool verbose_;

    bool verboseToFile_;

    int writePrecision_;

    word dirName_;

    mutable int rank_;

    mutable OFstream* sPtr;

    int printEvery_;

    bool sampleAll_;

    bool probeDebug_;

    bool includePosition_;

    const labelList particleIDsToSample_;

    mutable wordList itemsToSample_;

    mutable List<OFstream*> sPtrList_;                      //Must delete elements in destructor

    mutable int itemCounter_;

    mutable int currItemId_;

    mutable int printCounter_;

    mutable bool printNow_;

    mutable std::vector< std::vector<double> > sProbes_;
    
    mutable std::vector< std::vector<double*> > vProbes_;   //Must delete elements in destructor
    
    mutable  std::vector<std::string> probeName_;
    
    mutable int probeIndex_;
    
    mutable int printOnlyAt_;
  
    
public:

    //- Runtime type information
    TypeName("particleProbe");


    // Constructors

        //- Construct from components
        particleProbe
        (
            const dictionary& dict,
            cfdemCloud& sm,
            const word& typeName,
            const char* logFileName
        );

    // Destructor

        ~particleProbe();

    // Member Functions
    void updateProbes(int index, Field<scalar> sValues, Field<vector> vValues) const;
    void initialize(const word& modelName, const word& logFileName) const;
    void setOutputFile(const word& logFileName) const;
    void writeHeader() const;
    void writeProbe(int index, Field<scalar> sValues, Field<vector> vValues) const;
    bool checkIDForPrint(int) const;
    void setCounter() const;
    void clearProbes() const;
    std::vector< std::vector<double*> >* getVprobe() { return &vProbes_; };
    std::vector< std::vector<double> >* getSprobe() {return &sProbes_; };
    
    
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
