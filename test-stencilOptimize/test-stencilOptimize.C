/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of the IsoAdvector source code library, which is an
    unofficial extension to OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    isoSurf

Description
    Uses isoCutter to create a volume fraction field from either a cylinder,
    a sphere or a plane.

Author
    Henning Scheufler

\*---------------------------------------------------------------------------*/
#include "fvCFD.H"

#include "Field.H"
#include "fvMesh.H"
#include "volFields.H"
#include "OFstream.H"
#include "meshTools.H"

#include "extendedCentredCellToCellStencil.H"
#include "centredCPCCellToCellStencilObject.H"

#include "Random.H"
#include "zoneDistribute.H"



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"



    runTime.cpuTimeIncrement();

    // const extendedCentredCellToCellStencil& addressing =
    // centredCPCCellToCellStencilObject::New
    // (
    //     mesh
    // );

    Info << "building stencil took " << runTime.cpuTimeIncrement() << " s" << endl;


    runTime.cpuTimeIncrement();

    const stencil::extendedCentredCellToCellStencil& myaddressing =
    stencil::centredCPCCellToCellStencilObject::New
    (
        mesh
    );

    Info << "building mystencil took " << runTime.cpuTimeIncrement() << " s" << endl;






    return 0;
}


// ************************************************************************* //
