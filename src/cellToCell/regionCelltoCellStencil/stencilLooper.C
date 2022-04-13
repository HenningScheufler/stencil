/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2016 OpenFOAM Foundation
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

#include "stencilLooper.H"
#include "mapDistribute.H"
#include "cellToCellStencil.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::stencilLooper<Type>::stencilLooper
(
    const mapDistribute& map,
    const labelListList& stencil,
    const GeometricField<Type, Foam::fvPatchField, Foam::volMesh>& fld
)
:
    extendedCellToCellStencil(fld.mesh()),
    map_(map),
    stencil_(stencil),
    size_(stencil.size()),
    flatfld_(map_.constructSize())
{
    // Insert my internal values
    forAll(fld, celli)
    {
        flatfld_[celli] = fld[celli];
    }
    // Insert my boundary values
    label nCells = fld.mesh().nCells();
    label ninternalFaces = fld.mesh().nInternalFaces();
    forAll(fld.boundaryField(), patchi)
    {
        const fvPatchField<Type>& pfld = fld.boundaryField()[patchi];

        label nCompact = pfld.patch().start() - ninternalFaces + nCells;

        forAll(pfld, i)
        {
            flatfld_[nCompact++] = pfld[i];
        }
    }

    // Do all swapping
    map_.distribute(flatfld_);


}

template<class Type>
Foam::stencilLooper<Type>::stencilLooper
(
    const polyMesh& mesh,
    const mapDistribute& map,
    const labelListList& stencil,
    const Field<Type>& fld
)
:
    extendedCellToCellStencil(mesh),
    map_(map),
    stencil_(stencil),
    size_(stencil.size()),
    flatfld_(fld.size())
{
    // Insert my internal values
    forAll(fld, celli)
    {
        flatfld_[celli] = fld[celli];
    }

    // Do all swapping
    map_.distribute(flatfld_);


}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// template<class Type>
// void Foam::stencilLooper<Type>::compact()
// {
//     boolList isInStencil(map().constructSize(), false);
//     const labelListList& stencil = this->stencil();

//     forAll(stencil, celli)
//     {
//         const labelList& stencilCells = stencil[celli];

//         forAll(stencilCells, i)
//         {
//             isInStencil[stencilCells[i]] = true;
//         }
//     }

//     tmapPtr_.ref().compact(isInStencil, Pstream::msgType());
// }

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
const Foam::stencilLoop<Type> Foam::stencilLooper<Type>::operator[](const label i) const
{
    return stencilLoop<Type>(stencil_[i],flatfld_);
}

// ************************************************************************* //

template class Foam::stencilLooper<Foam::scalar>;
//template class Foam::stencilLooper<Foam::vector>;
