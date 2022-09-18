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

#ifndef HashTableSH_C
#define HashTableSH_C

#include "HashTableSH.H"
//#include "List.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

 
Foam::HashTableSH::HashTableSH(const label size)
:
    HashTableSHCore(),
    nElmts_(0),
    tableSize_(HashTableSHCore::canonicalSize(size)),
    table_(NULL)
{
    if (tableSize_)
    {
        table_ = new hashedEntry*[tableSize_];

        for (label hashIdx = 0; hashIdx < tableSize_; hashIdx++)
        {
            table_[hashIdx] = 0;
        }
    }
}


 
Foam::HashTableSH::HashTableSH(const HashTableSH& ht)
:
    HashTableSHCore(),
    nElmts_(0),
    tableSize_(ht.tableSize_),
    table_(NULL)
{
    if (tableSize_)
    {
        table_ = new hashedEntry*[tableSize_];

        for (label hashIdx = 0; hashIdx < tableSize_; hashIdx++)
        {
            table_[hashIdx] = 0;
        }

        for (const_iterator iter = ht.cbegin(); iter != ht.cend(); ++iter)
        {
            insert(iter.key(), *iter);
        }
    }
}

 /*
Foam::HashTableSH::HashTableSH
(
    const Xfer<HashTableSH>& ht
)
:
    HashTableSHCore(),
    nElmts_(0),
    tableSize_(0),
    table_(NULL)
{
    transfer(ht());
}
*/

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

 
Foam::HashTableSH::~HashTableSH()
{
    if (table_)
    {
        clear();
        delete[] table_;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

 
bool Foam::HashTableSH::found(const long int& key) const
{
    if (nElmts_)
    {
        const label hashIdx = hashKeyIndex(key);

        for (hashedEntry* ep = table_[hashIdx]; ep; ep = ep->next_)
        {
            if (key == ep->key_)
            {
                return true;
            }
        }
    }

    /* #ifdef FULLDEBUG
    if (debug)
    {
        InfoInFunction << "Entry " << key << " not found in hash table\n";
    }
    #endif */

    return false;
}


 
typename Foam::HashTableSH::iterator
Foam::HashTableSH::find
(
    const long int& key
)
{
    if (nElmts_)
    {
        const label hashIdx = hashKeyIndex(key);

        for (hashedEntry* ep = table_[hashIdx]; ep; ep = ep->next_)
        {
            if (key == ep->key_)
            {
                return iterator(this, ep, hashIdx);
            }
        }
    }

    /* #ifdef FULLDEBUG
    if (debug)
    {
        InfoInFunction << "Entry " << key << " not found in hash table\n";
    }
    #endif */

    return iterator();
}


 
typename Foam::HashTableSH::const_iterator
Foam::HashTableSH::find
(
    const long int& key
) const
{
    if (nElmts_)
    {
        const label hashIdx = hashKeyIndex(key);

        for (hashedEntry* ep = table_[hashIdx]; ep; ep = ep->next_)
        {
            if (key == ep->key_)
            {
                return const_iterator(this, ep, hashIdx);
            }
        }
    }

    /* #ifdef FULLDEBUG
    if (debug)
    {
        InfoInFunction << "Entry " << key << " not found in hash table\n";
    }
    #endif */

    return const_iterator();
}


 
 
bool Foam::HashTableSH::set
(
    const long int& key,
    const int& newEntry,
    const bool protect
)
{
    if (!tableSize_)
    {
        resize(2);
    }

    const label hashIdx = hashKeyIndex(key);

    hashedEntry* existing = 0;
    hashedEntry* prev = 0;

    for (hashedEntry* ep = table_[hashIdx]; ep; ep = ep->next_)
    {
        if (key == ep->key_)
        {
            existing = ep;
            break;
        }
        prev = ep;
    }

    // Not found, insert it at the head
    if (!existing)
    {
        table_[hashIdx] = new hashedEntry(key, table_[hashIdx], newEntry);
        nElmts_++;

        if (double(nElmts_)/tableSize_ > 0.8 && tableSize_ < maxTableSize)
        {
            /* #ifdef FULLDEBUG
            if (debug)
            {
                InfoInFunction << "Doubling table size\n";
            }
            #endif */

            resize(2*tableSize_);
        }
    }
    else if (protect)
    {
        // Found - but protected from overwriting
        // this corresponds to the STL 'insert' convention
        /* #ifdef FULLDEBUG
        if (debug)
        {
            InfoInFunction
                << "Cannot insert " << key << " already in hash table\n";
        }
        #endif */
        return false;
    }
    else
    {
        // Found - overwrite existing entry
        // this corresponds to the Perl convention
        hashedEntry* ep = new hashedEntry(key, existing->next_, newEntry);

        // Replace existing element - within list or insert at the head
        if (prev)
        {
            prev->next_ = ep;
        }
        else
        {
            table_[hashIdx] = ep;
        }

        delete existing;
    }

    return true;
}

// Shenghui 20191112		 
		bool Foam::HashTableSH::set
		(
			const long int& key,
			const int& newEntry,
			const bool protect,
			int& ZN
		)
		{		
		
		if (!tableSize_)
			{
				resize(2);
			}

			const label hashIdx = hashKeyIndex(key);

			hashedEntry* existing = 0;
			hashedEntry* prev = 0;

			for (hashedEntry* ep = table_[hashIdx]; ep; ep = ep->next_)
			{
				if (key == ep->key_)
				{
					existing = ep; 
					ZN = iterator(this, ep, hashIdx)();
					break;
				}
				prev = ep;
			}
			
			// Not found, insert it at the head
			if (!existing)
			{
				table_[hashIdx] = new hashedEntry(key, table_[hashIdx], newEntry);
				nElmts_++;
                //hashedEntry* ep = table_[hashIdx];
				
				if (double(nElmts_)/tableSize_> 0.8 && tableSize_ < maxTableSize)
				{
					/* #ifdef FULLDEBUG
					if (debug)
					{
						InfoInFunction << "Doubling table size\n";
					}
					#endif */

					this->resize(2*tableSize_);
				}
				
				//ZN = iterator(this, ep, hashIdx)();
			}
			else if (protect)
			{
				// Found - but protected from overwriting
				// this corresponds to the STL 'insert' convention
				/* #ifdef FULLDEBUG
				if (debug)
				{
					InfoInFunction
						<< "Cannot insert " << key << " already in hash table\n";
				}
				#endif */
				return false;
			}
			else
			{
				// Found - overwrite existing entry
				// this corresponds to the Perl convention
				hashedEntry* ep = new hashedEntry(key, existing->next_, newEntry);

				// Replace existing element - within list or insert at the head
				if (prev)
				{
					prev->next_ = ep;
				}
				else
				{
					table_[hashIdx] = ep;
				}

				delete existing;
			}

			return true;
		}		
//

 
bool Foam::HashTableSH::iteratorBase::erase()
{
    // Note: entryPtr_ is NULL for end(), so this catches that too
    if (entryPtr_)
    {
        // Search element before entryPtr_
        hashedEntry* prev = 0;

        for
        (
            hashedEntry* ep = hashTable_->table_[hashIndex_];
            ep;
            ep = ep->next_
        )
        {
            if (ep == entryPtr_)
            {
                break;
            }
            prev = ep;
        }

        if (prev)
        {
            // Has an element before entryPtr - reposition to there
            prev->next_ = entryPtr_->next_;
            delete entryPtr_;
            entryPtr_ = prev;
        }
        else
        {
            // entryPtr was first element on SLList
            hashTable_->table_[hashIndex_] = entryPtr_->next_;
            delete entryPtr_;

            // Assign any non-NULL pointer value so it doesn't look
            // like end()/cend()
            entryPtr_ = reinterpret_cast<hashedEntry*>(this);

            // Mark with special hashIndex value to signal it has been rewound.
            // The next increment will bring it back to the present location.
            //
            // From the current position 'curPos', we wish to continue at
            // prevPos='curPos-1', which we mark as markPos='-curPos-1'.
            // The negative lets us notice it is special, the extra '-1'
            // is needed to avoid ambiguity for position '0'.
            // To retrieve prevPos, we would later use '-(markPos+1) - 1'
            hashIndex_ = -hashIndex_ - 1;
        }

        hashTable_->nElmts_--;

        return true;
    }
    else
    {
        return false;
    }
}


 
bool Foam::HashTableSH::erase(const iterator& iter)
{
    // NOTE: We use (const iterator&) here, but manipulate its contents anyhow.
    // The parameter should be (iterator&), but then the compiler doesn't find
    // it correctly and tries to call as (iterator) instead.
    //
    // Adjust iterator after erase
    return const_cast<iterator&>(iter).erase();
}


 
bool Foam::HashTableSH::erase(const long int& key)
{
    return erase(find(key));
}


 


Foam::label Foam::HashTableSH::erase
(
    const HashTableSH& rhs
)
{
    label count = 0;

    // Remove rhs keys from this table - terminates early if possible
    // Could optimize depending on which hash is smaller ...
    for (iterator iter = begin(); iter != end(); ++iter)
    {
        if (rhs.found(iter.key()) && erase(iter))
        {
            count++;
        }
    }

    return count;
}


 
void Foam::HashTableSH::resize(const label sz)
{
    label newSize = HashTableSHCore::canonicalSize(sz);

    if (newSize == tableSize_)
    {
        /* #ifdef FULLDEBUG
          if (debug)
          {
            InfoInFunction << "New table size == old table size\n";
          }
        #endif */ 

        return;
    }

    HashTableSH* tmpTable = new HashTableSH(newSize);

    for (const_iterator iter = cbegin(); iter != cend(); ++iter)
    {
        tmpTable->insert(iter.key(), *iter);
    }

    label oldSize = tableSize_;
    tableSize_ = tmpTable->tableSize_;
    tmpTable->tableSize_ = oldSize;

    hashedEntry** oldTable = table_;
    table_ = tmpTable->table_;
    tmpTable->table_ = oldTable;

    delete tmpTable;
}


 
void Foam::HashTableSH::clear()
{
    if (nElmts_)
    {
        for (label hashIdx = 0; hashIdx < tableSize_; hashIdx++)
        {
            if (table_[hashIdx])
            {
                hashedEntry* ep = table_[hashIdx];
                while (hashedEntry* next = ep->next_)
                {
                    delete ep;
                    ep = next;
                }
                delete ep;
                table_[hashIdx] = 0;
            }
        }
        nElmts_ = 0;
    }
}


 
void Foam::HashTableSH::clearStorage()
{
    clear();
    resize(0);
}


 
void Foam::HashTableSH::shrink()
{
    const label newSize = HashTableSHCore::canonicalSize(nElmts_);

    if (newSize < tableSize_)
    {
        // Avoid having the table disappear on us
        resize(newSize ? newSize : 2);
    }
}


 
void Foam::HashTableSH::transfer(HashTableSH& ht)
{
    // As per the Destructor
    if (table_)
    {
        clear();
        delete[] table_;
    }

    tableSize_ = ht.tableSize_;
    ht.tableSize_ = 0;

    table_ = ht.table_;
    ht.table_ = NULL;

    nElmts_ = ht.nElmts_;
    ht.nElmts_ = 0;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

 
void Foam::HashTableSH::operator=
(
    const HashTableSH& rhs
)
{
    // Check for assignment to self
    if (this == &rhs)
    {
        FatalErrorInFunction
            << "attempted assignment to self"
            << abort(FatalError);
    }

    // Could be zero-sized from a previous transfer()
    if (!tableSize_)
    {
        resize(rhs.tableSize_);
    }
    else
    {
        clear();
    }

    for (const_iterator iter = rhs.cbegin(); iter != rhs.cend(); ++iter)
    {
        insert(iter.key(), *iter);
    }
}


 
bool Foam::HashTableSH::operator==
(
    const HashTableSH& rhs
) const
{
    // Sizes (number of keys) must match
    if (size() != rhs.size())
    {
        return false;
    }

    for (const_iterator iter = rhs.cbegin(); iter != rhs.cend(); ++iter)
    {
        const_iterator fnd = find(iter.key());

        if (fnd == cend() || fnd() != iter())
        {
            return false;
        }
    }

    return true;
}


 
bool Foam::HashTableSH::operator!=
(
    const HashTableSH& rhs
) const
{
    return !(operator==(rhs));
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

#include "HashTableSHIO.C"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
