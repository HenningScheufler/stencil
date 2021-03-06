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

#include "CPCCellToCellStencilUniqueListInsert.H"
#include "syncTools.H"
#include "dummyTransform.H"
#include "profiling.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::CPCCellToCellStencilUniqueListInsert::calcPointBoundaryData
(
    const boolList& isValidBFace,
    const labelList& boundaryPoints,
    Map<labelList>& neiGlobal
) const
{
    addProfiling(stencil, "CPCCellToCellStencilUniqueListInsert::calcPointBoundaryData");
    neiGlobal.resize(2*boundaryPoints.size());

    labelHashSet pointGlobals;

    forAll(boundaryPoints, i)
    {
        label pointi = boundaryPoints[i];

        neiGlobal.insert
        (
            pointi,
            calcFaceCells
            (
                isValidBFace,
                mesh().pointFaces()[pointi],
                pointGlobals
            )
        );
    }

    syncTools::syncPointMap
    (
        mesh(),
        neiGlobal,
        ListOps::unionEqOp(),
        Foam::dummyTransform()      // dummy transformation
    );
}


void Foam::CPCCellToCellStencilUniqueListInsert::calcCellStencil
(
    labelListList& globalCellCells
) const
{
    addProfiling(stencil, "CPCCellToCellStencilUniqueListInsert::calcCellStencil");
    // Calculate points on coupled patches
    labelList boundaryPoints(allCoupledFacesPatch()().meshPoints());


    // Mark boundary faces to be included in stencil (i.e. not coupled or empty)
    boolList isValidBFace;
    validBoundaryFaces(isValidBFace);


    // Swap pointCells for coupled points
    Map<labelList> neiGlobal;
    calcPointBoundaryData
    (
        isValidBFace,
        boundaryPoints,
        neiGlobal
    );

    globalCellCells.setSize(mesh().nCells());

    // Do coupled points first

    forAll(boundaryPoints, i)
    {
        label pointi = boundaryPoints[i];

        const labelList& pGlobals = neiGlobal[pointi];

        // Distribute to all pointCells
        const labelList& pCells = mesh().pointCells(pointi);

        forAll(pCells, j)
        {
            label celli = pCells[j];

            // Insert pGlobals into globalCellCells
            merge
            (
                globalNumbering().toGlobal(celli),
                pGlobals,
                globalCellCells[celli]
            );
        }
    }
    neiGlobal.clear();

    // Do remaining points cells
    labelHashSet pointGlobals;
    labelListList pGlobals(mesh().nPoints());
    DynamicList<label> faceCells(1000);

    for (label pointi = 0; pointi < mesh().nPoints(); pointi++)
    {
        calcFaceCells
        (
            isValidBFace,
            mesh().pointFaces()[pointi],
            faceCells // reset in function
        );

        pGlobals[pointi] = faceCells;
    }

    // oversized
    DynamicList<label> neiCells(1000);
    const labelListList& cellPoints = mesh().cellPoints();
    forAll(globalCellCells,celli)
    {
        neiCells.clear();
        const labelList& cPoints = cellPoints[celli];
        uniqueMerge(globalCellCells[celli],neiCells);
        for (const label pointi: cPoints)
        {
            uniqueMerge(pGlobals[pointi],neiCells);

        }
        moveIndexToStart(globalNumbering().toGlobal(celli),neiCells);

        globalCellCells[celli] = neiCells;

    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::CPCCellToCellStencilUniqueListInsert::CPCCellToCellStencilUniqueListInsert(const polyMesh& mesh)
:
    cellToCellStencil(mesh)
{
    // Calculate per cell the (point) connected cells (in global numbering)
    Info << " CPCCellToCellStencilUniqueListInsert " << endl;
    labelListList globalCellCells;
    calcCellStencil(*this);
}


// ************************************************************************* //
