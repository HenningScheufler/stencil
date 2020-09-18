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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::CPCCellToCellStencilUniqueListInsert::calcPointBoundaryData
(
    const boolList& isValidBFace,
    const labelList& boundaryPoints,
    Map<labelList>& neiGlobal
) const
{
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
        pGlobals[pointi] = calcFaceCells
        (
            isValidBFace,
            mesh().pointFaces()[pointi],
            pointGlobals
        );

        // pGlobals[pointi] = faceCells;
    }

    // oversized
    DynamicList<label> neiCells(1000);
    const labelListList& cellPoints = mesh().cellPoints();
    forAll(globalCellCells,celli)
    {
        neiCells.clear();
        // Info << "celli celli celli celli celli " << celli << endl;
        const labelList& cPoints = cellPoints[celli];
        neiCells.append(globalCellCells[celli]); // assume unique and sorted
        DynamicList<label>::iterator it_low;
        for (const label pointi: cPoints)
        {
            // neiCells.append(pGlobals[pointi]);
            // Info << "neiCells " << neiCells << endl;
            // Info << "pGlobals[pointi] " << pGlobals[pointi] << endl;
            label count = 0;
            for (const auto val: pGlobals[pointi])
            {
                // Info << "val " << val << "count " << count << endl;
                count++;
                if (neiCells.size() == 0)
                {
                    // Info << "no size" << endl;
                    neiCells.append(val);
                    continue; // break loop early
                }
                // Info << "neiCells " << neiCells << endl;
                it_low = std::lower_bound(neiCells.begin(), neiCells.end(), val);
                if (it_low == neiCells.end())
                {
                    // Info << "it low equals end" << endl;
                    neiCells.append(val);
                    // Info << "changed neiCells " << neiCells << endl;
                    continue; // break loop early
                }
                if (*it_low == val)
                {
                    // Info << "identical value " <<  val << endl;
                    // neiCells.append(val);
                    continue; // break loop early
                }

                // insert
                // Info << "insert " << endl;
                neiCells.setSize(neiCells.size()+1);
                // Info << "low " << std::distance(neiCells.begin(), it_low) << " val_low " << *it_low << endl;

                // it_low++;

                label start = std::distance(neiCells.begin(), it_low);
                // label tmp = *it_low;
                // *it_low = val;
                label tmp = val;
                // Info << "inserting val " << val << " " << " to neiCells " << neiCells  << endl;
                // Info << "tmp " << tmp << endl;
                // Info << "*it_low " << *it_low << endl;
                for (label i= start;i < neiCells.size();i++)
                {
                    std::swap(tmp,neiCells[i]);
                    // neiCells[i] = tmp;
                    // tmp = neiCells[i+1];
                }
                // Info << "inserting to neiCells " << neiCells << endl;
                // while(it_low != neiCells.end())
                // {
                    // Info << "it_low " << *it_low << " tmp " << tmp <<  endl;
                    // *i = tmp;
                    // tmp = *(i+1);
                    // *it_low = tmp;
                // }


                // neiCells.setSize(neiCells.size()+1);

                // for (DynamicList<label>::iterator i = it_low;i != neiCells.end();++i)
                // {
                //     Info << "iterator " << *i << endl;
                // }
                // Info << "low " << std::distance(neiCells.begin(), it_low) << " val_low " << *it_low << endl;

                // Info << "it pos " << it_low-neiCells.begin() << endl;
                // Info << "neiCells.end() " << neiCells.end() << endl;


                // neiCells.insert(low,val)
            }
        }

        globalCellCells[celli] = neiCells;

    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::CPCCellToCellStencilUniqueListInsert::CPCCellToCellStencilUniqueListInsert(const polyMesh& mesh)
:
    mycellToCellStencil(mesh)
{
    // Calculate per cell the (point) connected cells (in global numbering)
    labelListList globalCellCells;
    calcCellStencil(*this);
}


// ************************************************************************* //
