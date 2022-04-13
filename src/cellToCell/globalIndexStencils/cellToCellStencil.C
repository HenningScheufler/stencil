/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2018 OpenCFD Ltd.
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

#include "cellToCellStencil.H"
#include "syncTools.H"
#include "SortableList.H"
#include "emptyPolyPatch.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::cellToCellStencil::merge
(
    const label global0,
    const label global1,
    const labelList& listA,
    labelList& listB
)
{
    sort(listB);

    // See if global0, global1 already present in listB
    label nGlobalInsert = 0;

    if (global0 != -1)
    {
        label index0 = findSortedIndex(listB, global0);
        if (index0 == -1)
        {
            nGlobalInsert++;
        }
    }

    if (global1 != -1)
    {
        label index1 = findSortedIndex(listB, global1);
        if (index1 == -1)
        {
            nGlobalInsert++;
        }
    }


    // For all in listA see if they are present
    label nInsert = 0;

    forAll(listA, i)
    {
        label elem = listA[i];

        if (elem != global0 && elem != global1)
        {
            if (findSortedIndex(listB, elem) == -1)
            {
                nInsert++;
            }
        }
    }

    // Extend B with nInsert and whether global0,global1 need to be inserted.
    labelList result(listB.size() + nGlobalInsert + nInsert);

    label resultI = 0;

    // Insert global0,1 first
    if (global0 != -1)
    {
        result[resultI++] = global0;
    }
    if (global1 != -1)
    {
        result[resultI++] = global1;
    }


    // Insert listB
    forAll(listB, i)
    {
        label elem = listB[i];

        if (elem != global0 && elem != global1)
        {
            result[resultI++] = elem;
        }
    }


    // Insert listA
    forAll(listA, i)
    {
        label elem = listA[i];

        if (elem != global0 && elem != global1)
        {
            if (findSortedIndex(listB, elem) == -1)
            {
                result[resultI++] = elem;
            }
        }
    }

    if (resultI != result.size())
    {
        FatalErrorInFunction
            << "problem" << abort(FatalError);
    }

    listB.transfer(result);
}


void Foam::cellToCellStencil::uniqueMerge
(
    const labelList& listA,
    DynamicList<label>& sortedList
)
{
    // DynamicList<label>::iterator it_low;
    for (const auto val: listA)
    {
        addToSortedList(val,sortedList);
    }
}


void Foam::cellToCellStencil::addToSortedList
(
    const label idx,
    DynamicList<label>& sortedList
)
{
    DynamicList<label>::iterator it_low;

    if (sortedList.size() == 0)
    {
        sortedList.append(idx);
        return; // break loop early
    }

    it_low = std::lower_bound(sortedList.begin(), sortedList.end(), idx);
    if (it_low == sortedList.end())
    {
        sortedList.append(idx);
        return; // break loop early
    }
    if (*it_low == idx)
    {
        return; // break loop early
    }

    // insert
    sortedList.setSize(sortedList.size()+1);

    // label start = std::distance(sortedList.begin(), it_low);
    label tmp = idx;
    // for (label i= start;i < sortedList.size();i++)
    for (; it_low < sortedList.end(); it_low++)
    {
        std::swap(tmp,*it_low);
    }

}


void Foam::cellToCellStencil::moveIndexToStart
(
    const label globalId,
    List<label>& sortedList
)
{
    if (sortedList.size() <= 1)
    {
        return;
    }

    List<label>::iterator elem_it;
    elem_it = std::find(sortedList.begin(), sortedList.end(), globalId);
    if (elem_it == std::end(sortedList))
    {
        FatalErrorInFunction
            << "index not found in sortedList: "
            << sortedList
            << exit(FatalError);
    }


    List<label>::iterator it = sortedList.begin();
    label tmp = *it;
    *it = globalId;
    ++it;
    for (; it <= elem_it; ++it)
    {
        std::swap(tmp,*it);
    }

}


void Foam::cellToCellStencil::merge
(
    const label globalI,
    const labelList& pGlobals,
    labelList& cCells
)
{
    labelHashSet set;
    for (const label celli : cCells)
    {
        if (celli != globalI)
        {
            set.insert(celli);
        }
    }

    for (const label celli : pGlobals)
    {
        if (celli != globalI)
        {
            set.insert(celli);
        }
    }

    cCells.setSize(set.size()+1);
    label n = 0;
    cCells[n++] = globalI;

    for (const label seti : set)
    {
        cCells[n++] = seti;
    }
}


void Foam::cellToCellStencil::validBoundaryFaces(boolList& isValidBFace) const
{
    const polyBoundaryMesh& patches = mesh().boundaryMesh();

    isValidBFace.setSize(mesh().nBoundaryFaces(), true);

    for (const polyPatch& pp : patches)
    {
        if (pp.coupled() || isA<emptyPolyPatch>(pp))
        {
            label bFacei = pp.start()-mesh().nInternalFaces();
            forAll(pp, i)
            {
                isValidBFace[bFacei++] = false;
            }
        }
    }
}


Foam::autoPtr<Foam::indirectPrimitivePatch>
Foam::cellToCellStencil::allCoupledFacesPatch() const
{
    const polyBoundaryMesh& patches = mesh().boundaryMesh();

    label nCoupled = 0;

    for (const polyPatch& pp : patches)
    {
        if (pp.coupled())
        {
            nCoupled += pp.size();
        }
    }
    labelList coupledFaces(nCoupled);
    nCoupled = 0;

    for (const polyPatch& pp : patches)
    {
        if (pp.coupled())
        {
            label facei = pp.start();

            forAll(pp, i)
            {
                coupledFaces[nCoupled++] = facei++;
            }
        }
    }

    return autoPtr<indirectPrimitivePatch>::New
    (
        IndirectList<face>
        (
            mesh().faces(),
            coupledFaces
        ),
        mesh().points()
    );
}


void Foam::cellToCellStencil::insertFaceCells
(
    const label exclude0,
    const label exclude1,
    const boolList& isValidBFace,
    const labelList& faceLabels,
    labelHashSet& globals
) const
{
    const labelList& own = mesh().faceOwner();
    const labelList& nei = mesh().faceNeighbour();

    forAll(faceLabels, i)
    {
        label facei = faceLabels[i];

        label globalOwn = globalNumbering().toGlobal(own[facei]);
        if (globalOwn != exclude0 && globalOwn != exclude1)
        {
            globals.insert(globalOwn);
        }

        if (mesh().isInternalFace(facei))
        {
            label globalNei = globalNumbering().toGlobal(nei[facei]);
            if (globalNei != exclude0 && globalNei != exclude1)
            {
                globals.insert(globalNei);
            }
        }
        else
        {
            label bFacei = facei-mesh().nInternalFaces();

            if (isValidBFace[bFacei])
            {
                label globalI = globalNumbering().toGlobal
                (
                    mesh().nCells()
                  + bFacei
                );

                if (globalI != exclude0 && globalI != exclude1)
                {
                    globals.insert(globalI);
                }
            }
        }
    }
}


Foam::labelList Foam::cellToCellStencil::calcFaceCells
(
    const boolList& isValidBFace,
    const labelList& faceLabels,
    labelHashSet& globals
) const
{
    globals.clear();

    insertFaceCells
    (
        -1,
        -1,
        isValidBFace,
        faceLabels,
        globals
    );

    return globals.toc();
}


void Foam::cellToCellStencil::calcFaceCells
(
    const boolList& isValidBFace,
    const labelList& faceLabels,
    DynamicList<label>& globals
) const
{
    globals.clear();

    const labelList& own = mesh().faceOwner();
    const labelList& nei = mesh().faceNeighbour();
    const globalIndex& gblNum = globalNumbering();

    forAll(faceLabels, i)
    {
        label facei = faceLabels[i];

        label globalOwn = gblNum.toGlobal(own[facei]);
        addToSortedList(globalOwn,globals);

        if (mesh().isInternalFace(facei))
        {
            label globalNei = gblNum.toGlobal(nei[facei]);
            addToSortedList(globalNei,globals);
        }
        else
        {
            label bFacei = facei-mesh().nInternalFaces();

            if (isValidBFace[bFacei])
            {
                label globalI = gblNum.toGlobal
                (
                    mesh().nCells()
                  + bFacei
                );

                addToSortedList(globalI,globals);

            }
        }
    }

}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cellToCellStencil::cellToCellStencil(const polyMesh& mesh)
:
    mesh_(mesh),
    globalNumbering_(mesh_.nCells()+mesh_.nBoundaryFaces())
{}


// ************************************************************************* //
