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

#include <chrono>



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    // setup

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
        dimensionedScalar("0",dimless,100),
        "calculated"
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


    // building stencil
    Info<< "building stencil" << endl;
    runTime.cpuTimeIncrement();

    const extendedCentredCellToCellStencil& addressing =
    centredCPCCellToCellStencilObject::New
    (
        mesh
    );

    Info << "building CPC stencil took " << runTime.cpuTimeIncrement() << " s" << endl;

    scalar startZoneDistStencil  = std::clock();
    runTime.cpuTimeIncrement();

    zoneDistribute&  exchangeFields_ = zoneDistribute::New(mesh);

    // build whole stencil
    exchangeFields_.setUpCommforZone(boolList(mesh.nCells(),true));

    scalar endZoneDistStencil  = std::clock();
    Info << "build stencil zoneDistribute took " << runTime.cpuTimeIncrement() << " s" << endl;
    // Info << "build stencil zoneDistribute took " << endZoneDistStencil-startZoneDistStencil << " µs" << endl;

    // Info<< "cellCellCell:" << endl;
    // writeStencilStats(addressing.stencil());

    // Collect stencil cell centres
    Info<< "Collecting data 100 times" << endl;
    runTime.cpuTimeIncrement();
    scalar startStdStencil  = std::clock();

    scalar countCollectData = 0;
    scalar countLooping = 0;


    scalar count = 0;

    for(label i = 0;i<100;i++)
    {
        countCollectData -=  std::clock();
        List<List<scalar>> stencilCellNums(mesh.nCells());

        addressing.collectData
        (
            cellNumbers,
            stencilCellNums
        );
        countCollectData +=  std::clock();


        count = 0;
        countLooping -=  std::clock();
        forAll(selected,celli)
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
        countLooping +=  std::clock();
    }

    scalar endStdStencil  = std::clock();

    Info << "standard stencil took to add all values" << runTime.cpuTimeIncrement() << " s" << endl;
    Info << "standard stencil took to add all values " << endStdStencil-startStdStencil << " µs" << endl;
    Info << "standard stencil took to add all values " << endStdStencil-startStdStencil << " µs" << endl;
    Info << "countCollectData " << countCollectData << " µs" << endl;
    Info << "countLooping " << countLooping << " µs" << endl;


    Info << "standard stencil count " << count << endl;
    Info << " " << endl;

    // zoneDistribute
    runTime.cpuTimeIncrement();

    scalar startZoneDist = std::clock();
    countCollectData = 0;
    countLooping = 0;
    for(label i = 0;i<100;i++)
    {
        countCollectData -=  std::clock();
        exchangeFields_.setUpCommforZone(selected,false);

        const labelListList& stencil = exchangeFields_.getStencil();


        Map<scalar> mapCellNum =
            exchangeFields_.getDatafromOtherProc(selected,cellNumbers);
        countCollectData +=  std::clock();

        countLooping -=  std::clock();
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

        countLooping +=  std::clock();

    }
    scalar endZoneDist  = std::clock();

    Info << "zoneDistribute took to add all values " << runTime.cpuTimeIncrement() << " s" << endl;
    Info << "zoneDistribute took to add all values " << endZoneDist-startZoneDist << " µs" << endl;
    Info << "countCollectData " << countCollectData << " µs" << endl;
    Info << "countLooping " << countLooping << " µs" << endl;

    Info << "zoneDistribute count " << count << endl;


    // best case scenario looping
    runTime.cpuTimeIncrement();


    List<List<scalar>> stencilCellNums(mesh.nCells());
    label nSelectedCells = 0;

    for (auto sel: selected)
    {
        if(sel)
        {
            nSelectedCells++;
        }
    }
    Info << "nSelectedCells " << nSelectedCells << endl;
    Info << "mesh.nCells() " << mesh.nCells() << endl;

    addressing.collectData
    (
        cellNumbers,
        stencilCellNums
    );

    List<List<scalar>> selectedStencilCells(nSelectedCells);
    forAll(selected,celli)
    {
        if(selected[celli])
        {
            selectedStencilCells.append(stencilCellNums[celli]);
        }
    }
    scalar startBestCase = std::clock();
    countCollectData = 0;
    countLooping = 0;
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();


    for(label i = 0;i<100;i++)
    {

        count = 0;
        countLooping -=  std::clock();
        for (const auto& stencilVal: selectedStencilCells)
        {
            for (const auto& val: stencilVal)
            {
                count += val;
            }
        }

        reduce(count,sumOp<scalar>());
        countLooping +=  std::clock();

    }
    scalar endBestCase  = std::clock();

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

    std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[µs]" << std::endl;
    std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::nanoseconds> (end - begin).count() << "[ns]" << std::endl;

    Info << "BestCase took to add all values " << runTime.cpuTimeIncrement() << " s" << endl;
    Info << "BestCase took to add all values " << endBestCase-startBestCase << " µs" << endl;
    Info << "BestCase countCollectData " << countCollectData << " µs" << endl;
    Info << "BestCase countLooping " << countLooping << " µs" << endl;

    Info << "BestCase count " << count << endl;

    // adding all values
    scalar startAddingAll = std::clock();
    runTime.cpuTimeIncrement();

    for(label i = 0;i<100;i++)
    {

        count = 0;
        countLooping -=  std::clock();
        for (const auto& stencilVal: stencilCellNums)
        {
            for (const auto& val: stencilVal)
            {
                count += val;
            }
        }

        reduce(count,sumOp<scalar>());
        countLooping +=  std::clock();

    }

    scalar endAddingAll = std::clock();
    Info << "Adding all took to add all values " << runTime.cpuTimeIncrement() << " s" << endl;
    Info << "Adding all values countLooping " << endAddingAll-startAddingAll  << " µs" << endl;


        // adding all values
    scalar startBitSet = std::clock();
    runTime.cpuTimeIncrement();
    labelList cellis = maskCells.toc();

    for(label i = 0;i<100;i++)
    {

        count = 0;
        countLooping -=  std::clock();
        // for (const auto& stencilVal: stencilCellNums)
        // {
        //     for (const auto& val: stencilVal)
        //     {
        //         count += val;
        //     }
        // }
        // forAll(selected,celli)
        for (const label& celli: cellis)
        {
            // if(selected[celli])
            // {
                for (const auto val: stencilCellNums[celli])
                {
                    count += val;
                }
            // }
        }

        reduce(count,sumOp<scalar>());
        countLooping +=  std::clock();

    }

    scalar endBitSet = std::clock();
    Info << "Adding bitset took to add all values " << runTime.cpuTimeIncrement() << " s" << endl;
    Info << "Adding bitSet values countLooping " << endBitSet-startBitSet  << " µs" << endl;



    profiling::writeNow();



    return 0;
}


// ************************************************************************* //
