/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
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

#include "error.H"

// * * * * * * * * * * * * * Private Member Classes * * * * * * * * * * * * //


inline Foam::HashTableSH::hashedEntry::hashedEntry
(
    const long int& key,
    hashedEntry* next,
    const int& obj
)
:
    key_(key),
    next_(next),
    obj_(obj)
{}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //


inline Foam::label
Foam::HashTableSH::hashKeyIndex(const long int& key) const
{
    // size is power of two - this is the modulus
    return string::hash()(key) & (tableSize_ - 1);
}


// * * * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * //


inline Foam::label Foam::HashTableSH::capacity() const
{
    return tableSize_;
}



inline Foam::label Foam::HashTableSH::size() const
{
    return nElmts_;
}



inline bool Foam::HashTableSH::empty() const
{
    return !nElmts_;
}



inline bool Foam::HashTableSH::insert
(
    const long int& key,
    const int& newEntry
)
{
    return this->set(key, newEntry, true);
}



inline bool Foam::HashTableSH::set
(
    const long int& key,
    const int& newEntry
)
{
    return this->set(key, newEntry, false);
}


/*
inline Foam::Xfer<Foam::HashTableSH>
Foam::HashTableSH::xfer()
{
    return xferMove(*this);
}
*/

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //


inline int& Foam::HashTableSH::operator[](const long int& key)
{
    iterator iter = this->find(key);

    if (iter == this->end())
    {
        FatalErrorInFunction
            << key << " not found in table.  Valid entries: "
          //  << toc()
            << exit(FatalError);
    }

    return *iter;
}



inline const int& Foam::HashTableSH::operator[](const long int& key) const
{
    const_iterator iter = this->find(key);

    if (iter == this->cend())
    {
        FatalErrorInFunction
            << key << " not found in table.  Valid entries: "
            //<< toc()
            << exit(FatalError);
    }

    return *iter;
}



inline int& Foam::HashTableSH::operator()(const long int& key)
{
    iterator iter = this->find(key);

    if (iter == this->end())
    {
        this->insert(key, int());
        return *find(key);
    }
    else
    {
        return *iter;
    }
}


// * * * * * * * * * * * * * * * iterator base * * * * * * * * * * * * * * * //


inline Foam::HashTableSH::iteratorBase::iteratorBase()
:
    hashTable_(0),
    entryPtr_(0),
    hashIndex_(0)
{}



inline Foam::HashTableSH::iteratorBase::iteratorBase
(
    const HashTableSH* hashTbl
)
:
    hashTable_(const_cast<HashTableSH*>(hashTbl)),
    entryPtr_(0),
    hashIndex_(0)
{
    if (hashTable_->nElmts_)
    {
        // find first non-NULL table entry
        while
        (
            !(entryPtr_ = hashTable_->table_[hashIndex_])
         && ++hashIndex_ < hashTable_->tableSize_
        )
        {}

        if (hashIndex_ >= hashTable_->tableSize_)
        {
            // make into an end iterator
            entryPtr_ = 0;
            hashIndex_ = 0;
        }
    }
}



inline Foam::HashTableSH::iteratorBase::iteratorBase
(
    const HashTableSH* hashTbl,
    const hashedEntry* elmt,
    const label hashIndex
)
:
    hashTable_(const_cast<HashTableSH*>(hashTbl)),
    entryPtr_(const_cast<hashedEntry*>(elmt)),
    hashIndex_(hashIndex)
{}



inline void
Foam::HashTableSH::iteratorBase::increment()
{
    // A negative index is a special value from erase
    if (hashIndex_ < 0)
    {
        // the markPos='-curPos-1', but we wish to continue at 'curPos-1'
        // thus use '-(markPos+1) -1'
        hashIndex_ = -(hashIndex_+1) - 1;
    }
    else if (entryPtr_)
    {
        if (entryPtr_->next_)
        {
            // Move to next element on the SLList
            entryPtr_ = entryPtr_->next_;
            return;
        }
    }
    // else
    // {
    //     // if we reach here (entryPtr_ is NULL) it is already at the end()
    //     // we should probably stop
    // }


    // Step to the next table entry
    while
    (
        ++hashIndex_ < hashTable_->tableSize_
     && !(entryPtr_ = hashTable_->table_[hashIndex_])
    )
    {}

    if (hashIndex_ >= hashTable_->tableSize_)
    {
        // make into an end iterator
        entryPtr_ = 0;
        hashIndex_ = 0;
    }
}



inline
const long int& Foam::HashTableSH::iteratorBase::key() const
{
    return entryPtr_->key_;
}



inline int&
Foam::HashTableSH::iteratorBase::object()
{
    return entryPtr_->obj_;
}



inline const int&
Foam::HashTableSH::iteratorBase::cobject() const
{
    return entryPtr_->obj_;
}



inline bool Foam::HashTableSH::iteratorBase::operator==
(
    const iteratorBase& iter
) const
{
    return entryPtr_ == iter.entryPtr_;
}



inline bool Foam::HashTableSH::iteratorBase::operator!=
(
    const iteratorBase& iter
) const
{
    return entryPtr_ != iter.entryPtr_;
}



inline bool Foam::HashTableSH::iteratorBase::operator==
(
    const iteratorEnd&
) const
{
    return !entryPtr_;
}



inline bool Foam::HashTableSH::iteratorBase::operator!=
(
    const iteratorEnd&
) const
{
    return entryPtr_;
}


// * * * * * * * * * * * * * * * * STL iterator  * * * * * * * * * * * * * * //


inline Foam::HashTableSH::iterator::iterator()
:
    iteratorBase()
{}



inline Foam::HashTableSH::iterator::iterator
(
    const iteratorEnd&
)
:
    iteratorBase()
{}



inline Foam::HashTableSH::iterator::iterator
(
    HashTableSH* hashTbl
)
:
    iteratorBase(hashTbl)
{}



inline Foam::HashTableSH::iterator::iterator
(
    HashTableSH* hashTbl,
    hashedEntry* elmt,
    const label hashIndex
)
:
    iteratorBase(hashTbl, elmt, hashIndex)
{}



inline int&
Foam::HashTableSH::iterator::operator*()
{
    return this->object();
}



inline int&
Foam::HashTableSH::iterator::operator()()
{
    return this->object();
}



inline const int&
Foam::HashTableSH::iterator::operator*() const
{
    return this->cobject();
}



inline const int&
Foam::HashTableSH::iterator::operator()() const
{
    return this->cobject();
}



inline
typename Foam::HashTableSH::iterator&
Foam::HashTableSH::iterator::operator++()
{
    this->increment();
    return *this;
}



inline typename Foam::HashTableSH::iterator
Foam::HashTableSH::iterator::operator++(int)
{
    iterator old = *this;
    this->increment();
    return old;
}



inline typename Foam::HashTableSH::iterator
Foam::HashTableSH::begin()
{
    return iterator(this);
}


// * * * * * * * * * * * * * * * STL const_iterator * * * * * * * * * * * * * //


inline Foam::HashTableSH::const_iterator::const_iterator()
:
    iteratorBase()
{}



inline Foam::HashTableSH::const_iterator::const_iterator
(
    const HashTableSH::iterator& iter
)
:
    iteratorBase(iter)
{}



inline Foam::HashTableSH::const_iterator::const_iterator
(
    const iteratorEnd&
)
:
    iteratorBase()
{}



inline Foam::HashTableSH::const_iterator::const_iterator
(
    const HashTableSH* hashTbl
)
:
    iteratorBase(hashTbl)
{}



inline Foam::HashTableSH::const_iterator::const_iterator
(
    const HashTableSH* hashTbl,
    const hashedEntry* elmt,
    const label hashIndex
)
:
    iteratorBase(hashTbl, elmt, hashIndex)
{}



inline const int&
Foam::HashTableSH::const_iterator::operator*() const
{
    return this->cobject();
}



inline const int&
Foam::HashTableSH::const_iterator::operator()() const
{
    return this->cobject();
}



inline
typename Foam::HashTableSH::const_iterator&
Foam::HashTableSH::const_iterator::operator++()
{
    this->increment();
    return *this;
}



inline typename Foam::HashTableSH::const_iterator
Foam::HashTableSH::const_iterator::operator++(int)
{
    const_iterator old = *this;
    this->increment();
    return old;
}



inline typename Foam::HashTableSH::const_iterator
Foam::HashTableSH::cbegin() const
{
    return const_iterator(this);
}



inline typename Foam::HashTableSH::const_iterator
Foam::HashTableSH::begin() const
{
    return this->cbegin();
}


// ************************************************************************* //
