/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

\*---------------------------------------------------------------------------*/

#include "maskedCPCStencil.H"
#include "centredCPCCellToCellStencilObject.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


void Foam::maskedCPCStencil::calcCellStencil
(
    labelListList& globalCellCells
)
{
    const extendedCentredCellToCellStencil& stdCPCstencil =
        centredCPCCellToCellStencilObject::New(mesh_);


    label fieldSize = stdCPCstencil.map().constructSize();

    labelList addressing(fieldSize,-1);

    // simpler without BCs
    bitSet inStencil(fieldSize,false);
    label nCells = mesh_.nCells();
    if (includeBoundaryCells_)
    {
        nCells += mesh_.nBoundaryFaces();
    }

    // is in compact number not global numbering
    const labelListList& stencil = stdCPCstencil.stencil();

    // the stencil is only defined at the cells not the interface
    // but interface values can be included in the stencil
    label nStencilCells = 0;
    for (const label celli:celliAddressing_)
    {
        if (celli >= mesh_.nCells())
        {
            ++nStencilCells;
            break;
        }
        for (const label nei:stencil[celli])
        {
            if (nei < nCells && mask_.test(celli)) // only the local cells
            {
                inStencil.set(nei);
            }
        }
    }

    globalNumbering_ = globalIndex(inStencil.count());
    label localSize = 0;
    for (const label idx:inStencil)
    {
        addressing[idx] = globalNumbering_.toGlobal(localSize);
        ++localSize;
    }

    stdCPCstencil.map().distribute(addressing);


    globalCellCells.setSize(nStencilCells);
    DynamicList<label> subSetNeiCells(1000);

    localSize = 0;
    for (label i=0;i < nStencilCells;i++)
    {
        const label celli = celliAddressing_[i];
        subSetNeiCells.clear();
        for (const label nei:stencil[celli])
        {
            const label newIdx = addressing[nei];
            if (newIdx != -1)
            {
                subSetNeiCells.append(newIdx);
            }
        }
        globalCellCells[localSize] = subSetNeiCells;
        ++localSize;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


Foam::maskedCPCStencil::maskedCPCStencil
(
    const fvMesh& mesh,
    const bitSet& mask,
    bool includeBoundaryCells
)
:
    mesh_(mesh),
    globalNumbering_(),
    map_(nullptr),
    mask_(mask),
    celliAddressing_(mask.toc()),
    includeBoundaryCells_(includeBoundaryCells)
{
    // Calculate per cell the (point) connected cells (in global numbering)
    // labelListList globalCellCells;
    calcCellStencil(*this);

    List<Map<label>> compactList;
    map_.reset
    (
        new mapDistribute
        (
            globalNumbering_,
            *this,
            compactList
        )
    );

    List<label> test= identity(map_().constructSize());

    map_().distribute(test);
}


// ************************************************************************* //
