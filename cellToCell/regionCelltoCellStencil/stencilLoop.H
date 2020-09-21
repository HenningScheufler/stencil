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

Class
    Foam::stencilLoop

Description

SourceFiles
    stencilLoop.C

\*---------------------------------------------------------------------------*/

#ifndef stencilLoop_H
#define stencilLoop_H

#include "labelList.H"
#include "volFieldsFwd.H"
#include "volFields.H"
#include "GeometricField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators
// template<class Type>
// class GeometricField<Type, Foam::fvPatchField, Foam::volMesh>;

// template<class Type>
// inline const Type& operator[](const label i);

/*---------------------------------------------------------------------------*\
              Class stencilLoop Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class stencilLoop
:
    public labelList
{
    // Private data

    const GeometricField<Type, Foam::fvPatchField, Foam::volMesh>& fld_;

    // Private Member Functions

        //- No copy construct
        stencilLoop
        (
            const stencilLoop&
        ) = delete;

        //- No copy assignment
        void operator=(const stencilLoop&) = delete;


public:

    // Forward declaration of access classes

        class const_iterator;

    // Constructors

        //- Construct from uncompacted cell stencil
        // explicit stencilLoop(const cellToCellStencil&);
        stencilLoop
        (
            const labelList& neiCells,
            const GeometricField<Type, Foam::fvPatchField, Foam::volMesh>& fld
        );


// Iteration

        //- A const_iterator for iterating across \a on values
        class const_iterator
        {
            // template<class Type>
            friend class stencilLoop; //<Type>;

            //- The parent being iterated
            const stencilLoop* stencil_;

            //- Global position of the current \a on bit
            label pos_;

            //- Default construct - an end iterator
            inline const_iterator() noexcept;

            //- Construct begin iterator
            inline const_iterator(const stencilLoop* loop);

            //- Construct iterator, starting at or beyond the given position
            inline const_iterator(const stencilLoop* loop, label pos);

        public:

            //- Return the current \a on position
            inline const Type& operator*() const noexcept;

            //- Move to the next \a on position
            inline const_iterator& operator++();

            inline bool operator==(const const_iterator& iter) const noexcept;
            inline bool operator!=(const const_iterator& iter) const noexcept;
        };


        //- Iterator set to the position of the first \a on bit
        inline const_iterator begin() const;

        //- Iterator set to the position of the first \a on bit
        inline const_iterator cbegin() const;

        //- Iterator set to the position of the first \a on bit that occurs
        //- at or beyond the given position
        inline const_iterator begin(label pos) const;

        //- Iterator set to the position of the first \a on bit that occurs
        //- at or beyond the given position
        inline const_iterator cbegin(label pos) const;

        //- Iterator beyond the end of the bitSet
        inline const_iterator end() const noexcept;

        //- Iterator beyond the end of the bitSet
        inline const_iterator cend() const noexcept;

    // Member Functions

        const labelList& indices()
        {
            return *this;
        }


        const Type& operator[](const label i) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#include "stencilLoopI.H"


#endif

// ************************************************************************* //
