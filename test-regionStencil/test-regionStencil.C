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

#include "stencilLoop.H"
#include "regionCelltoCellStencil.H"
#include "CPCCellToCellStencil.H"
#include "extendedCelltoCellStencilLooper.H"
#include "extendedCentredCellToCellStencil.H"
#include "centredCPCCellToCellStencilObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

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
        cellNumbers[celli] = globalNumbering.toGlobal(celli) + 100;
    }

    // labelList test = identity(mesh.nCells());
    // stencilLoop<scalar> loop(test,cellNumbers);

    // // test operator
    // forAll(loop,i)
    // {
    //     Info << "loop operator [] " << loop[i] << endl;
    //     Info << "indices " << loop.indices()[i] << endl;
    // }

    // // test iterators
    // for (auto& a:loop)
    // {
    //     Info<< "non const loop " << a << endl;
    // }

    // for (const auto& a:loop)
    // {
    //     Info<< "const loop " << a << endl;
    // }

    runTime.cpuTimeIncrement();

    CPCCellToCellStencil addressing(mesh);

    Info << "building uniqueListsInsert stencil took " << runTime.cpuTimeIncrement() << " s" << endl;
    boolList selectedCells(mesh.nCells(),false);

    if (Pstream::master())
    {
        selectedCells = true;
    }

    Map<scalar> mapTest;
    mapTest.insert(1,20.0);

    // Info << " addressing " << addressing << endl;

    const extendedCentredCellToCellStencil& stdCPCstencil =
        centredCPCCellToCellStencilObject::New(mesh);

            // tmp<mapDistribute>& tmapPtr_,
            // tmp<labelListList>& tstencil_,
            // const GeometricField<Type, Foam::fvPatchField, Foam::volMesh>& fld
    // tmp<mapDistribute> tmapPtr (new mapDistribute(stdCPCstencil.map()));
    // tmp<labelListList> tstencil (new labelListList(stdCPCstencil.stencil()));

    extendedCelltoCellStencilLooper<scalar> stencilLooper
    (
        stdCPCstencil.map(),
        stdCPCstencil.stencil(),
        cellNumbers
    );

        // 1. Construct cell data in compact addressing
    List<scalar> flatFld(stencilLooper.map().constructSize(), Zero);

    // Insert my internal values
    forAll(cellNumbers, celli)
    {
        flatFld[celli] = cellNumbers[celli];
    }
    // Insert my boundary values
    forAll(cellNumbers.boundaryField(), patchi)
    {
        const fvPatchField<scalar>& pfld = cellNumbers.boundaryField()[patchi];

        label nCompact =
            pfld.patch().start()
           -cellNumbers.mesh().nInternalFaces()
           +cellNumbers.mesh().nCells();

        forAll(pfld, i)
        {
            flatFld[nCompact++] = pfld[i];
        }
    }

    // Do all swapping
    stencilLooper.map().distribute(flatFld);

    // Info << "flatFld " << flatFld << endl;

    // Info << "stencilLoop flatFld " << stencilLoop.flatField() << endl;


    List<List<scalar>> stencilValues;
    stencilLooper.collectData
    (
        stencilValues
    );

    // Info << " stencilValues " << stencilValues << endl;

    // labelList stencilASd = stencilLoop[0];
    // Field<scalar> testFld(flatFld);
    // labelList test = identity(mesh.nCells())

    // stencilLoop<scalar> loop(test,cellNumbers);
    // stencilLoop<scalar> looptest(test,testFld);

    forAll(stencilLooper,celli)
    {
        stencilLoop<scalar> loop = stencilLooper[celli];

        // forAll(loop,i)
        // {
        //     Info << "loop operator [] " << loop[i] << endl;
        //     Info << "indices " << loop.indices()[i] << endl;
        // }

        for (const auto& a:loop)
        {
            Info<< "const loop " << a << endl;
        }
    }




    return 0;
}


// ************************************************************************* //