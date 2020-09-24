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

#include "CPCCellToCellStencilCellBased.H"
#include "CPCCellToCellStencil.H"
#include "CPCCellToCellStencilUniqueList.H"
#include "CPCCellToCellStencilUniqueListInsert.H"
#include "CPCCellToCellStencilCellBasedBitSet.H"



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
void checkAddressing
(
    const labelListList& stencil1,
    const labelListList& stencil2
)
{
    bool correct = true;
    forAll(stencil1,celli)
    {
        bool cellCorrect = true;
        if (stencil1[celli][0] != stencil2[celli][0])
        {
            correct = false;
            cellCorrect = false;
        }
        for (const label val:stencil1[celli])
        {
            if(!stencil2[celli].found(val))
            {
                correct = false;
                cellCorrect = false;
            }
        }
        if (!cellCorrect)
        {
            // Pout << "do not match celli " << celli << endl;
            // Pout << "stencil2 " << stencil2[celli]  << endl;
            // Pout << "stencil1                 " << stencil1[celli]  << endl;
        }
    }
    if (correct)
    {
        Pout << "stencil is correct" << endl;
    }
    else
    {
        Pout << "stencil are not match" << endl;
    }

}
int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"



    runTime.cpuTimeIncrement();

    CPCCellToCellStencil addressing(mesh);

    Info << "building stencil took " << runTime.cpuTimeIncrement() << " s" << endl;
    // Info << "addressing " << addressing << endl;


    runTime.cpuTimeIncrement();

    CPCCellToCellStencilCellBased cellBasedaddressing(mesh);

    Info << "building  cellbased merge stencil took " << runTime.cpuTimeIncrement() << " s" << endl;
    // Info << "cellBasedaddressing " << cellBasedaddressing << endl;


    runTime.cpuTimeIncrement();

    CPCCellToCellStencilUniqueList uniqueListaddressing(mesh);

    Info << "building uniqueLists stencil took " << runTime.cpuTimeIncrement() << " s" << endl;
    // Info << "uniqueListaddressing " << uniqueListaddressing << endl;



    runTime.cpuTimeIncrement();

    CPCCellToCellStencilUniqueListInsert uniqueInsertListaddressing(mesh);

    Info << "building uniqueListsInsert stencil took " << runTime.cpuTimeIncrement() << " s" << endl;
    // Info << "uniqueInsertListaddressing " << uniqueInsertListaddressing << endl;

    runTime.cpuTimeIncrement();

    CPCCellToCellStencilCellBasedBitSet bitSetAddressing(mesh);

    Info << "building bitSetAddressing stencil took " << runTime.cpuTimeIncrement() << " s" << endl;

    checkAddressing(addressing,cellBasedaddressing);
    checkAddressing(addressing,uniqueListaddressing);
    checkAddressing(addressing,uniqueInsertListaddressing);



    profiling::writeNow();

    OFstream profData("profiling_buildStencil.dat");
    Foam::profiling::print(profData);

    return 0;
}


// ************************************************************************* //
