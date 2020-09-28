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


// Iterators
template<class Type>
Foam::stencilLooper<Type> Foam::subSetCPCStencil::stencilValues
(
    const GeometricField<Type, Foam::fvPatchField, Foam::volMesh>& fld
) const
{

    const labelList& send = this->sendIndices();
    Field<scalar> flatField(send.size(),0);

    label nCells = mesh_.nCells();

    forAll(send,i) //const label celli:maskedStencil.celliAddressing())
    {
        label idx = send[i];
        if (idx == -1)
        {
            continue;
        }

        if (i < nCells)
        {
            flatField[idx] = fld[i];
        }
        else
        {
            const label faceI = i + mesh_.nInternalFaces() - mesh_.nCells();
            const polyBoundaryMesh& pbm = mesh_.boundaryMesh();
            // Boundary face. Find out which face of which patch
            const label patchI = pbm.whichPatch(faceI);
            if (patchI < 0 || patchI >= pbm.size())
            {
                    FatalErrorInFunction
                    << "Cannot find patch for face " << faceI
                    << abort(FatalError);
            }
            const polyPatch& pp = pbm[patchI];

            if (fld.boundaryField()[patchI].size())
            {
                const label patchFaceI = pp.whichFace(faceI);
                flatField[idx] = fld.boundaryField()[patchI][patchFaceI];
            }

        }

    }

    return stencilLooper<Type>
    (
        mesh_,
        map(),
        *this,
        flatField
    );
}

// ************************************************************************* //
