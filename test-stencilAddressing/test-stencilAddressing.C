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
#include "profiling.H"
#include "extendedCelltoCellStencilLooper.H"

#include "subSetCPCStencil.H"
#include "maskedCPCStencil.H"

#include "maskedCPCStencilCompact.H"



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void loopStdStencil
(
    const boolList& selected,
    const volScalarField& cellNum,
    const extendedCentredCellToCellStencil& addressing,
    const label n,
    word profDesc
)
{
    const fvMesh& mesh = cellNum.mesh();
    // stops if it gets out of scope;
    profilingTrigger loopStencil(profDesc);
    scalar count = 0;

        profilingTrigger getStencilValues(profDesc + "buildStencil");
        List<List<scalar>> stencilCellNums(mesh.nCells());

        addressing.collectData
        (
            cellNum,
            stencilCellNums
        );
        getStencilValues.stop();


    for(label i = 0;i<n;i++)
    {
        count = 0;
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
    }
    Info << "std stencil count is " << count << endl;
}


void loopZoneDist
(
    const boolList& selected,
    const volScalarField& cellNum,
    zoneDistribute& exchangeFields_,
    const label n,
    word profDesc
)
{
    const fvMesh& mesh = cellNum.mesh();
    profilingTrigger loopStencil(profDesc);
    scalar count = 0;


        profilingTrigger getStencilValues(profDesc + "buildStencil");
        exchangeFields_.setUpCommforZone(selected,false);

        const labelListList& stencil = exchangeFields_.getStencil();



        Map<scalar> mapCellNum =
            exchangeFields_.getDatafromOtherProc(selected,cellNum);

        getStencilValues.stop();

    for(label i = 0;i<n;i++)
    {
        count = 0;
        forAll(selected,celli)
        {
            if (selected[celli])
            {
                for (const label gblIdx : stencil[celli])
                {
                    count += exchangeFields_.getValue(cellNum, mapCellNum, gblIdx);
                }
            }
        }
        reduce(count,sumOp<scalar>());

    }
    Info << "zoneDistribute count is " << count << endl;
}

void loopStencilLooper
(
    const boolList& selected,
    const volScalarField& cellNum,
    const extendedCentredCellToCellStencil& addressing,
    const label n,
    word profDesc
)
{
    // stops if it gets out of scope;
    profilingTrigger loopStencil(profDesc);
    scalar count = 0;


        profilingTrigger getStencilValues(profDesc + "buildStencil");
        extendedCelltoCellStencilLooper<scalar> stencilLooper
        (
            addressing.map(),
            addressing.stencil(),
            cellNum
        );
        getStencilValues.stop();

    for(label i = 0;i<n;i++)
    {
        count = 0;
        forAll(selected,celli)
        {
            if(selected[celli])
            {
                // stencilLoop<scalar> loop = stencilLooper[celli];
                for (const auto& val:stencilLooper[celli])
                {
                    count += val;
                }

            }
        }

        reduce(count,sumOp<scalar>());
    }
    Info << "stencil loop count is " << count << endl;
}

void subSetStencil
(
    const boolList& selected,
    const volScalarField& cellNum,
    const extendedCentredCellToCellStencil& addressing,
    const label n,
    word profDesc
)
{
    // stops if it gets out of scope;
    const fvMesh& mesh = cellNum.mesh();
    bitSet cells(selected);
    // cells.set(selected);
    labelList neededCells  = cells.toc();
    profilingTrigger loopStencil(profDesc);
    scalar count = 0;


    profilingTrigger getStencilValues(profDesc + "buildStencil");
    subSetCPCStencil subsetStencil(cellNum.mesh(),neededCells,true);
        // extendedCelltoCellStencilLooper<scalar> stencilLooper
        // (
        //     addressing.map(),
        //     addressing.stencil(),
        //     cellNum
        // );
    getStencilValues.stop();

    for(label i = 0;i<n;i++)
    {

        const labelList& cellAddress = subsetStencil.subSetCells();
        const labelList& sendIndices = subsetStencil.sendIndices();
        Field<scalar> cellNumField(mesh.nCells()+mesh.nBoundaryFaces(),0);

        label nCells = cellNum.mesh().nCells();

        // Info << " sendIndices size " << sendIndices.size() << " constructSize.size() " << subsetStencil.map().constructSize() << endl;
        // Info << " cellNumField size " << cellNumField.size() <<  endl;

        // Info << " sendIndices  " << sendIndices << endl;

        forAll(sendIndices,i) //const label celli:maskedStencil.celliAddressing())
        {
            label idx = sendIndices[i];
            if (idx == -1)
            {
                continue;
            }

            if (i < nCells)
            {
                cellNumField[idx] = cellNum[i];
            }
            else
            {
                const label faceI = i + mesh.nInternalFaces() - mesh.nCells();

                const polyBoundaryMesh& pbm = mesh.boundaryMesh();

                // Boundary face. Find out which face of which patch
                const label patchI = pbm.whichPatch(faceI);

                if (patchI < 0 || patchI >= pbm.size())
                {
                        FatalErrorInFunction
                        << "Cannot find patch for face " << faceI
                        << abort(FatalError);
                }
                const polyPatch& pp = pbm[patchI];
                // Info << "pp.name " << pp.name() << endl;
                if (isA<emptyPolyPatch>(pp))
                {
                    cellNumField[idx] = 0;
                }
                else
                {
                    if (cellNum.boundaryField()[patchI].size())
                    {
                        const label patchFaceI = pp.whichFace(faceI);
                        cellNumField[idx] = cellNum.boundaryField()[patchI][patchFaceI];
                        // label asd = cellNum.boundaryField()[patchI][patchFaceI];
                    }
                }


            }

            // stencilLoop<scalar> loop = stencilLooper[celli];
        }

        subsetStencil.map().distribute(cellNumField);


        count = 0;
        for(const auto& neiCelli:subsetStencil)
        {
            for (const label idx:neiCelli)
            {
                count += cellNumField[idx];
            }
        }
        // Info << "not counting " << endl;

        reduce(count,sumOp<scalar>());
    }
    Info << "subsetStencil loop count is " << count << endl;
}


void maskMaskedStencil
(
    const boolList& selected,
    const volScalarField& cellNum,
    const extendedCentredCellToCellStencil& addressing,
    const label n,
    word profDesc
)
{
    // stops if it gets out of scope;
    const fvMesh& mesh = cellNum.mesh();
    // bitSet cells2(selected.size(),true);
    bitSet cells(cellNum.mesh().nCells()+cellNum.mesh().nBoundaryFaces(),false);
    // Info << "cellNum" << cellNum << endl;
    forAll(selected,celli)
    {
        if (selected[celli])
        {
            cells.set(celli);
        }
    }


    // labelList neededCells  = cells2.toc();
    profilingTrigger loopStencil(profDesc);
    scalar count = 0;

    for(label i = 0;i<n;i++)
    {

        profilingTrigger getStencilValues(profDesc + "buildStencil");
        maskedCPCStencilCompact maskedStencil(cellNum.mesh(),cells,true);

        // const labelList& cellAddress = maskedStencil.subSetCells();
        const labelList& sendIndices = maskedStencil.sendIndices();
        Field<scalar> cellNumField(sendIndices.size(),0);

        label nCells = cellNum.mesh().nCells();

        forAll(sendIndices,i) //const label celli:maskedStencil.celliAddressing())
        {
            label idx = sendIndices[i];
            if (idx == -1)
            {
                continue;
            }

            if (i < nCells)
            {
                cellNumField[idx] = cellNum[i];
            }
            else
            {
                const label faceI = i + mesh.nInternalFaces() - mesh.nCells();

                const polyBoundaryMesh& pbm = mesh.boundaryMesh();

                // Boundary face. Find out which face of which patch
                const label patchI = pbm.whichPatch(faceI);

                if (patchI < 0 || patchI >= pbm.size())
                {
                        FatalErrorInFunction
                        << "Cannot find patch for face " << faceI
                        << abort(FatalError);
                }
                const polyPatch& pp = pbm[patchI];
                // Info << "pp.name " << pp.name() << endl;
                if (isA<emptyPolyPatch>(pp))
                {
                    cellNumField[idx] = 0;
                }
                else
                {
                    if (cellNum.boundaryField()[patchI].size())
                    {
                        const label patchFaceI = pp.whichFace(faceI);
                        cellNumField[idx] = cellNum.boundaryField()[patchI][patchFaceI];
                        // label asd = cellNum.boundaryField()[patchI][patchFaceI];
                    }
                }
            }


            // stencilLoop<scalar> loop = stencilLooper[celli];
        }

        maskedStencil.map().distribute(cellNumField);



        count = 0;
        for(const auto& neiCelli:maskedStencil)
        {
            for (const label idx:neiCelli)
            {
                count += cellNumField[idx];
            }
        }

        getStencilValues.stop();



        reduce(count,sumOp<scalar>());
    }
    Info << "maskMaskedStencil loop count is " << count << endl;
}

void maskedStdStencil
(
    const boolList& selected,
    const volScalarField& cellNum,
    const extendedCentredCellToCellStencil& addressing,
    const label n,
    word profDesc
)
{
    const fvMesh& mesh = cellNum.mesh();

    bitSet cells(cellNum.mesh().nCells()+cellNum.mesh().nBoundaryFaces(),false);
    // Info << "cellNum" << cellNum << endl;
    forAll(selected,celli)
    {
        if (selected[celli])
        {
            cells.set(celli);
        }
    }

    List<bool> maskedCells = cells.values();
    // stops if it gets out of scope;
    profilingTrigger loopStencil(profDesc);
    scalar count = 0;

    profilingTrigger getStencilValues(profDesc + "buildStencil");
    List<List<scalar>> stencilCellNums(mesh.nCells());

    addressing.collectData
    (
        cellNum,
        stencilCellNums
    );

    addressing.map().distribute(maskedCells);
    getStencilValues.stop();


    for(label i = 0;i<n;i++)
    {
        count = 0;
        forAll(selected,celli)
        {
            if(selected[celli])
            {
                forAll(addressing.stencil()[celli],i)
                {
                    const label idx = addressing.stencil()[celli][i];
                    if (maskedCells[idx])
                    {
                        count += stencilCellNums[celli][i];
                    }
                }
            }
        }

        reduce(count,sumOp<scalar>());
    }
    Info << "masked std stencil count is " << count << endl;
}


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
    profilingTrigger profBuildStencil("profBuildStencil");


    const extendedCentredCellToCellStencil& addressing =
    centredCPCCellToCellStencilObject::New
    (
        mesh
    );

    zoneDistribute&  exchangeFields_ = zoneDistribute::New(mesh);

    // build whole stencil
    exchangeFields_.setUpCommforZone(boolList(mesh.nCells(),true));

    profBuildStencil.stop();

    // Collect stencil cell centres
    Info<< "Collecting data 100 times" << endl;

    loopStdStencil(selected,cellNumbers,addressing,1,word("partialStdStencil"));

    loopZoneDist(selected,cellNumbers,exchangeFields_,1,word("partialZoneDist"));

    loopStencilLooper(selected,cellNumbers,addressing,1,word("partialStencilLooper"));

    subSetStencil(selected,cellNumbers,addressing,1,word("partialsubSetStencil"));

    // masked loop

    selectCells.setSize(mesh.nCells()*0.3);
    forAll(selectCells, i)
    {
        selectCells[i] = rndGen.position<label>(0,mesh.nCells()-1);
    }

    maskCells = false;

    maskCells.set(selectCells);

    selected = maskCells.values();
    Info << "masked"  << endl;

    maskedStdStencil(selected,cellNumbers,addressing,1,word("maskedStdStencil"));

    maskMaskedStencil(selected,cellNumbers,addressing,1,word("maskedmaskedStencil"));

    selected = true;
    Info << "full"  << endl;

    loopStdStencil(selected,cellNumbers,addressing,1,word("fullStdStencil"));

    loopZoneDist(selected,cellNumbers,exchangeFields_,1,word("fullZoneDist"));

    loopStencilLooper(selected,cellNumbers,addressing,1,word("fullStencilLooper"));

    subSetStencil(selected,cellNumbers,addressing,1,word("fullsubSetStencil"));

    //maskMaskedStencil(selected,cellNumbers,addressing,1,word("fullmaskedStencil"));



    OFstream profData("profiling_stencilAddressing.dat");
    Foam::profiling::print(profData);



    return 0;
}


// ************************************************************************* //
