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
#include "profiling.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


void Foam::maskedCPCStencil::calcCellStencil
(
    labelListList& globalCellCells
)
{
    const extendedCentredCellToCellStencil& stdCPCstencil =
        centredCPCCellToCellStencilObject::New(mesh_);



    const labelListList& stencil = stdCPCstencil.stencil();

    List<Map<label>> compactList;
    map_.reset
    (
        new mapDistribute
        (
            stdCPCstencil.map()
        )
    );




    boolList isInStencil(map().constructSize(), false);

    for (const label idx: mask_)
    {
        isInStencil[idx] = true;
    }

    stdCPCstencil.map().distribute(isInStencil);

    labelList oldToNewSub;
    labelList oldToNewConstruct;
    map_.ref().compact
    (
        isInStencil,
        mesh_.nCells()+mesh_.nBoundaryFaces(),      // maximum index of subMap
        sendMapIndices_,
        oldToNewConstruct,
        UPstream::msgType()
    );


    DynamicList<label> neiCells(100);
    globalCellCells.setSize(celliAddressing_.size());
    label count= 0;
    for (const label celli:celliAddressing_)
    {
        neiCells.clear();
        for (const label nei:stencil[celli])
        {
            label newIdx = oldToNewConstruct[nei];
            if (newIdx != -1 && isInStencil[nei])
            {
                neiCells.append(newIdx);
            }
        }
        globalCellCells[count] = neiCells;
        ++count;
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
    map_(nullptr),
    mask_(mask),
    celliAddressing_(mask.toc()),
    includeBoundaryCells_(includeBoundaryCells)
{
    if (mask.size() != mesh.nCells()+mesh.nBoundaryFaces())
    {
        FatalError  << "the size mask is to be expected to be " << nl
                    << "mesh.nCells()+mesh.nBoundaryFaces()" << nl
                    << "if the boundary values or not needed set" << nl
                    << "mesh.nBoundaryFaces() to false" << nl
                    << abort(FatalError);
    }

    // Calculate per cell the (point) connected cells (in global numbering)
    // labelListList globalCellCells;
    labelListList& stencil = *this;
    calcCellStencil(stencil); // map sendMapIndices_ is set here
}

// ************************************************************************* //
