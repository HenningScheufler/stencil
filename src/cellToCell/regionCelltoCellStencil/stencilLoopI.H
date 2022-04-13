/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2019 OpenCFD Ltd.
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


// * * * * * * * * * * * * * * * * Iterators * * * * * * * * * * * * * * * * //

template<class Type>
inline Foam::stencilLoop<Type>::const_iterator::const_iterator() noexcept
:
    stencil_(nullptr),
    pos_(0)
{}

template<class Type>
inline Foam::stencilLoop<Type>::const_iterator::const_iterator(const stencilLoop* loop)
:
    stencil_(loop),
    pos_(0)
{}

template<class Type>
inline Foam::stencilLoop<Type>::const_iterator::const_iterator
(
    const stencilLoop* loop,
    label pos
)
:
    stencil_(loop),
    pos_(pos)
{}

template<class Type>
inline const Type& Foam::stencilLoop<Type>::const_iterator::operator*() const noexcept
{
    return stencil_->operator[](pos_);
}

template<class Type>
inline typename Foam::stencilLoop<Type>::const_iterator& Foam::stencilLoop<Type>::const_iterator::operator++()
{
    pos_++;
    return *this;
}

template<class Type>
inline bool Foam::stencilLoop<Type>::const_iterator::operator==
(
    const const_iterator& iter
) const noexcept
{
    return (iter.pos_ == pos_);
}

template<class Type>
inline bool Foam::stencilLoop<Type>::const_iterator::operator!=
(
    const const_iterator& iter
) const noexcept
{
    return (iter.pos_ != pos_);
}

template<class Type>
inline typename Foam::stencilLoop<Type>::const_iterator Foam::stencilLoop<Type>::begin() const
{
    return const_iterator(this);
}

template<class Type>
inline typename Foam::stencilLoop<Type>::const_iterator Foam::stencilLoop<Type>::cbegin() const
{
    return const_iterator(this);
}

template<class Type>
inline typename Foam::stencilLoop<Type>::const_iterator Foam::stencilLoop<Type>::begin(label pos) const
{
    return const_iterator(this, pos);
}

template<class Type>
inline typename Foam::stencilLoop<Type>::const_iterator Foam::stencilLoop<Type>::cbegin(label pos) const
{
    return const_iterator(this, pos);
}

template<class Type>
inline typename Foam::stencilLoop<Type>::const_iterator Foam::stencilLoop<Type>::end() const noexcept
{
    return const_iterator(this,neiIdx_.size());
}

template<class Type>
inline typename Foam::stencilLoop<Type>::const_iterator Foam::stencilLoop<Type>::cend() const noexcept
{
    return const_iterator(this,neiIdx_.size());
}
