/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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
    Foam::centredListBasedCPCCellToCellStencilObject

Description

SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef centredListBasedCPCCellToCellStencilObject_H
#define centredListBasedCPCCellToCellStencilObject_H

#include "extendedCentredCellToCellStencil.H"
#include "CPCCellToCellStencilUniqueListInsert.H"
#include "MeshObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class centredListBasedCPCCellToCellStencilObject Declaration
\*---------------------------------------------------------------------------*/

class centredListBasedCPCCellToCellStencilObject
:
    public MeshObject
    <
        fvMesh,
        TopologicalMeshObject,
        centredListBasedCPCCellToCellStencilObject
    >,
    public extendedCentredCellToCellStencil
{

public:

    TypeName("centredListBasedCPCCellToCellStencil");

    // Constructors

        //- Construct from uncompacted cell stencil
        explicit centredListBasedCPCCellToCellStencilObject
        (
            const fvMesh& mesh
        )
        :
            MeshObject
            <
                fvMesh,
                Foam::TopologicalMeshObject,
                centredListBasedCPCCellToCellStencilObject
            >(mesh),
            extendedCentredCellToCellStencil(CPCCellToCellStencilUniqueListInsert(mesh))
        {}


    //- Destructor
    virtual ~centredListBasedCPCCellToCellStencilObject() = default;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
