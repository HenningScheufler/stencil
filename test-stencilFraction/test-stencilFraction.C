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

    bitSet maskCells(mesh.nCells(),false);

    Random rndGen(123456);// identical to default seed

    labelList selectCells(label(mesh.nCells()*0.01));
    forAll(selectCells, i)
    {
        selectCells[i] = rndGen.position<label>(0,mesh.nCells()-1);
    }

    maskCells.set(selectCells);

    boolList selected = maskCells.values();


    runTime.cpuTimeIncrement();

    const extendedCentredCellToCellStencil& addressing =
    centredCPCCellToCellStencilObject::New
    (
        mesh
    );

    Info << "building stencil took " << runTime.cpuTimeIncrement() << " s" << endl;

    // Info<< "cellCellCell:" << endl;
    // writeStencilStats(addressing.stencil());

    // Collect stencil cell centres
    runTime.cpuTimeIncrement();
    scalar startStdStencil  = std::clock();

    List<List<scalar>> stencilCellNums(mesh.nCells());
    addressing.collectData
    (
        cellNumbers,
        stencilCellNums
    );

    scalar count = 0;
    forAll(stencilCellNums,celli)
    {
        if(selected[celli])
        {
            for (const auto val: stencilCellNums[celli])
            {
                count += val;
            }
        }
    }

    reduce(count,sumOp<scalar>());

    scalar endStdStencil  = std::clock();

    Info << "standard stencil took " << runTime.cpuTimeIncrement() << " s" << endl;
    Info << "standard stencil took " << endStdStencil-startStdStencil << " ms" << endl;
    Info << "standard stencil count " << count << endl;

    // zoneDistribute

    scalar startZoneDistStencil  = std::clock();
    runTime.cpuTimeIncrement();

    zoneDistribute&  exchangeFields_ = zoneDistribute::New(mesh);

    // build whole stencil
    exchangeFields_.setUpCommforZone(boolList(mesh.nCells(),true));

    scalar endZoneDistStencil  = std::clock();
    Info << "build stencil zoneDistribute took " << runTime.cpuTimeIncrement() << " s" << endl;
    Info << "build stencil zoneDistribute took " << endZoneDistStencil-startZoneDistStencil << " ms" << endl;

    scalar startZoneDist = std::clock();

    exchangeFields_.setUpCommforZone(selected);

    const labelListList& stencil = exchangeFields_.getStencil();

    Map<scalar> mapCellNum =
        exchangeFields_.getDatafromOtherProc(selected,cellNumbers);

    count = 0;
    forAll(selected,celli)
    {
        if (selected[celli])
        {
            for (const label gblIdx : stencil[celli])
            {
                count += exchangeFields_.getValue(cellNumbers, mapCellNum, gblIdx);
            }

        }
    }
    reduce(count,sumOp<scalar>());
    scalar endZoneDist  = std::clock();

    Info << "zoneDistribute took " << endZoneDist-startZoneDist << " ms" << endl;
    Info << "zoneDistribute count " << count << endl;








    return 0;
}


// ************************************************************************* //
