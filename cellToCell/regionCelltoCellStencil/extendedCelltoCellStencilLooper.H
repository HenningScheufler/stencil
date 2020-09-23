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
    Foam::extendedCelltoCellStencilLooper

Description

SourceFiles
    extendedCelltoCellStencilLooper.C

\*---------------------------------------------------------------------------*/

#ifndef extendedCelltoCellStencilLooper_H
#define extendedCelltoCellStencilLooper_H

#include "extendedCellToCellStencil.H"
#include "extendedCellToFaceStencil.H"
#include "stencilLoop.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class cellToCellStencil;

/*---------------------------------------------------------------------------*\
              Class extendedCelltoCellStencilLooper Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class extendedCelltoCellStencilLooper
:
    public extendedCellToCellStencil
{
    // Private data

        //- Swap map for getting neighbouring data
        const mapDistribute& map_;

        //- Per cell the stencil.
        const labelListList& stencil_;

        const label size_;

        const GeometricField<Type, Foam::fvPatchField, Foam::volMesh>& fld_;

        Field<Type> flatfld_;


    // Private Member Functions

        //- No copy construct
        extendedCelltoCellStencilLooper
        (
            const extendedCelltoCellStencilLooper&
        ) = delete;

        //- No copy assignment
        void operator=(const extendedCelltoCellStencilLooper&) = delete;


public:

    // Constructors

        //- Construct from uncompacted cell stencil
        explicit extendedCelltoCellStencilLooper
        (
            const mapDistribute& map,
            const labelListList& stencil,
            const GeometricField<Type, Foam::fvPatchField, Foam::volMesh>& fld
        );


    // Member Functions

        // Iteration
        label size()
        {
            return size_;
        }

        //- Return reference to the parallel distribution map
        const mapDistribute& map() const
        {
            return map_;
        }

        //- Return reference to the stencil
        const labelListList& stencil() const
        {
            return stencil_;
        }

        const GeometricField<Type, Foam::fvPatchField, Foam::volMesh>& geoField()
        {
            return fld_;
        }

        const Field<Type>& flatField()
        {
            return flatfld_;
        }

        // - After removing elements from the stencil adapt the schedule (map).
        // void compact();

        //- Use map to get the data into stencil order
        void collectData
        (
            List<List<Type>>& stencilFld
        ) const
        {
            extendedCellToFaceStencil::collectData
            (
                map(),
                stencil(),
                fld_,
                stencilFld
            );
        }

        //- Use map to get the data into stencil order
        // List<stencilLoop<Type>> collectData
        // (
        //     const GeometricField<Type, fvPatchField, volMesh>& fld
        // ) const
        // {

        // }

        //- Sum vol field contributions to create cell values
        template<class WeightType>
        tmp
        <
            GeometricField
            <
                typename outerProduct<WeightType, Type>::type,
                fvPatchField,
                volMesh
            >
        > weightedSum
        (
            const GeometricField<Type, fvPatchField, volMesh>& fld,
            const List<List<WeightType>>& stencilWeights
        ) const
        {
            return extendedCellToCellStencil::weightedSum
            (
                map(),
                stencil(),
                fld,
                stencilWeights
            );
        }

        const stencilLoop<Type> operator[](const label i) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
