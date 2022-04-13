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

#include "stencilLoop.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::stencilLoop<Type>::stencilLoop
(
    const labelList& neiCells,
    const Field<Type>& fld
)
:
    neiIdx_(neiCells),
    size_(neiCells.size()),
    fld_(fld)
{

}

template<class Type>
Foam::stencilLoop<Type>::stencilLoop
(
    const stencilLoop& loop
)
:
    neiIdx_(loop.neiIdx_),
    size_(loop.size_),
    fld_(loop.fld_)
{

}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
const Type& Foam::stencilLoop<Type>::operator[](const label i) const
{
    return fld_[neiIdx_[i]];
}
// ************************************************************************* //


template class Foam::stencilLoop<Foam::scalar>;
// template class Foam::stencilLoop<Foam::vector>;