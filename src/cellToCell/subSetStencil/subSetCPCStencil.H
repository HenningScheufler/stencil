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

Class
    Foam::subSetCPCStencil

Description

SourceFiles
    subSetCPCStencil.C

\*---------------------------------------------------------------------------*/

#ifndef subSetCPCStencil_H
#define subSetCPCStencil_H

#include "boolList.H"
#include "HashSet.H"
#include "Map.H"
#include "fvMesh.H"
#include "globalIndex.H"
#include "mapDistribute.H"

#include "stencilLooper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class subSetCPCStencil Declaration
\*---------------------------------------------------------------------------*/

class subSetCPCStencil
:
    public labelListList
{
    // Private data

        const fvMesh& mesh_;

        autoPtr<mapDistribute> map_;

        const labelList& subSetCells_;

        labelList sendMapIndices_;

        const bool includeBoundaryCells_;



        //- Calculates per cell the neighbour data
        //  (= cell or boundary in global numbering).
        //  First element is always cell itself!
        void calcCellStencil(labelListList& globalCellCells);


public:

    // Constructors

        //- Construct from all cells and boundary faces
        explicit subSetCPCStencil
        (
            const fvMesh& mesh,
            const labelList& subSetCells,
            bool includeBoundaryCells
        );


        const mapDistribute& map() const
        {
            return map_();
        }

        // const bitSet& mask()
        // {
        //     return mask_;
        // }
        const labelList& sendIndices() const
        {
            return sendMapIndices_;
        }

        const labelList& subSetCells() const
        {
            return subSetCells_;
        }

        // Iterators
        template<class Type>
        stencilLooper<Type> stencilValues
        (
            const GeometricField<Type, Foam::fvPatchField, Foam::volMesh>& fld
        ) const;

        template<class Type>
        stencilLooper<Type> stencilValues
        (
            const Field<Type>& fld
        ) const
        {
            const labelListList& stencil = *this;
            return stencilLooper<Type>
            (
                mesh_,
                this->map(),
                stencil,
                fld
            );
        }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#include "subSetCPCStencilI.C"

#endif

// ************************************************************************* //
