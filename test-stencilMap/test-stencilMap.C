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

#include "CPCCellToCellStencil.H"
#include "CPCCellToCellStencilCellBased.H"
#include "CPCCellToCellStencilUniqueList.H"
#include "CPCCellToCellStencilUniqueListInsert.H"



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
template<class Type>
void collectData
(
    const mapDistribute& map,
    const labelListList& stencil,
    const GeometricField<Type, fvPatchField, volMesh>& fld,
    List<List<Type>>& stencilFld
)
{
    Info << "stencil " << stencil << endl;
    // 1. Construct cell data in compact addressing
    List<Type> flatFld(map.constructSize(), Zero);

    // Insert my internal values
    forAll(fld, celli)
    {
        flatFld[celli] = fld[celli];
    }
    // Insert my boundary values
    forAll(fld.boundaryField(), patchi)
    {
        const fvPatchField<Type>& pfld = fld.boundaryField()[patchi];

        label nCompact =
            pfld.patch().start()
           -fld.mesh().nInternalFaces()
           +fld.mesh().nCells();

        forAll(pfld, i)
        {
            flatFld[nCompact++] = pfld[i];
        }
    }
    Info << "flatFld " << flatFld << endl;

    // Do all swapping
    map.distribute(flatFld);

    Info << "flatFld " << flatFld << endl;

    // 2. Pull to stencil
    stencilFld.setSize(stencil.size());

    forAll(stencil, facei)
    {
        const labelList& compactCells = stencil[facei];
        Info << "compactCells " << compactCells << endl;

        stencilFld[facei].setSize(compactCells.size());

        forAll(compactCells, i)
        {
            stencilFld[facei][i] = flatFld[compactCells[i]];
        }
    }
}

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

    volScalarField cellNumbers
    (
        IOobject
        (
            "cellNumbers",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("0",dimless,0),
        "zeroGradient"
    );
    globalIndex globalNumbering(mesh.nCells());

    forAll(cellNumbers,celli)
    {
        cellNumbers[celli] = globalNumbering.toGlobal(celli);
    }



    runTime.cpuTimeIncrement();

    CPCCellToCellStencil addressing(mesh);

    Info << "building stencil took " << runTime.cpuTimeIncrement() << " s" << endl;
    // Info << "addressing " << addressing << endl;

    // runTime.cpuTimeIncrement();

    // CPCCellToCellStencilUniqueListInsert uniqueInsertListaddressing(mesh);

    // Info << "building uniqueListsInsert stencil took " << runTime.cpuTimeIncrement() << " s" << endl;
    // Info << "uniqueInsertListaddressing " << uniqueInsertListaddressing << endl;


    // checkAddressing(addressing,cellBasedaddressing);
    // checkAddressing(addressing,uniqueListaddressing);
    // checkAddressing(addressing,uniqueInsertListaddressing);
//
    // Pout << "addressing.globalNumbering().offSets() " << addressing.globalNumbering().offsets() << endl;
    // Pout << "addressing.globalNumbering().sizes() " << addressing.globalNumbering().sizes() << endl;
    Info << "addressing " << addressing << endl;

   autoPtr<mapDistribute> mapPtr_;

    // Calculate distribute map (also renumbers elements in stencil)
    List<Map<label>> compactMap(Pstream::nProcs());
    mapPtr_.reset
    (
        new mapDistribute
        (
            addressing.globalNumbering(),
            addressing,
            compactMap
        )
    );
    Info << "addressing " << addressing << endl;
    Info << "compactMap " << compactMap << endl;
    // List<List<vector>> test(mesh.nCells());
    List<List<scalar>> test(mesh.nCells());

    // Pout << "test " << test.size() << endl;

    collectData
    (
        mapPtr_(),
        addressing,
        cellNumbers,
        // mesh.C(),
        test
    );





    return 0;
}


// ************************************************************************* //