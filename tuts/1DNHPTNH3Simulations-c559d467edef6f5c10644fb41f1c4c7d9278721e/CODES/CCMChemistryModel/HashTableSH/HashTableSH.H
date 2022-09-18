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

Class
    Foam::HashTableSH

Description
    An STL-conforming hash table.

Note
    Hashing index collisions are handled via chaining using a singly-linked
    list with the colliding entry being added to the head of the linked
    list. Thus copying the hash table (or indeed even resizing it) will
    often result in a different hash order. Use a sorted table-of-contents
    when the hash order is important.

SourceFiles
    HashTableSHI.H
    HashTableSH.C
    HashTableSHIO.C

\*---------------------------------------------------------------------------*/

#ifndef HashTableSH_H
#define HashTableSH_H

#include "label.H"
//#include "uLabel.H"
#include "word.H"
#include "className.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

//template<int> class List;
//template<int> class UList;
//template<long int> class List;
//template<long int> class UList;
class HashTableSH;
template<class T, class Key, class Hash> class HashPtrTable;


Istream& operator>>(Istream&, HashTableSH&);


Ostream& operator<<(Ostream&, const HashTableSH&);


/*---------------------------------------------------------------------------*\
                         Class HashTableSHCore Declaration
\*---------------------------------------------------------------------------*/

//- Template-invariant bits for HashTableSH
struct HashTableSHCore
{
    //- Return a canonical (power-of-two) size
    static label canonicalSize(const label);

    //- Maximum allowable table size
    static const label maxTableSize;

    //- Construct null
    HashTableSHCore()
    {}

    //- Define template name and debug
   // ClassName("HashTableSH");

    //- A zero-sized end iterator
    struct iteratorEnd
    {
        //- Construct null
        iteratorEnd()
        {}
    };

    //- iteratorEnd set to beyond the end of any HashTableSH
    inline static iteratorEnd cend()
    {
        return iteratorEnd();
    }

    //- iteratorEnd set to beyond the end of any HashTableSH
    inline static iteratorEnd end()
    {
        return iteratorEnd();
    }
};


/*---------------------------------------------------------------------------*\
                           Class HashTableSH Declaration
\*---------------------------------------------------------------------------*/

//template<class int, class long int=word, class Hash=string::hash>
class HashTableSH
:
    public HashTableSHCore
{
    // Private data type for table entries

        //- Structure to hold a hashed entry with SLList for collisions
        struct hashedEntry
        {
            //- The lookup key
            long int key_;

            //- Pointer to next hashedEntry in sub-list
            hashedEntry* next_;

            //- The data object
            int obj_;

            //- Construct from key, next pointer and object
            inline hashedEntry(const long int&, hashedEntry* next, const int&);


        private:
            //- Disallow default bitwise copy construct
            hashedEntry(const hashedEntry&);

            //- Disallow default bitwise assignment
            void operator=(const hashedEntry&);
        };


    // Private data: size of table, the table and current number of elements

        //- The current number of elements in table
        label nElmts_;

        //- Number of primary entries allocated in table
        label tableSize_;

        //- The table of primary entries
        hashedEntry** table_;


    // Private Member Functions

        //- Return a canonical (power-of-two) size
        static label canonicalSize(const label);

        //- Return the hash index of the long int within the current table size.
        //  No checks for zero-sized tables.
        inline label hashKeyIndex(const long int&) const;

        //- Assign a new hashedEntry to a possibly already existing key
        bool set(const long int&, const int& newElmt, bool protect);
		


public:

    // Forward declaration of iterators

        class iteratorBase;
        class iterator;
        class const_iterator;

        //- Declare friendship with the HashPtrTable class
        template<class T2, class Key2, class Hash2>
        friend class HashPtrTable;

        //- Declare friendship with the iteratorBase
        friend class iteratorBase;

        //- Declare friendship with the iterator
        friend class iterator;

        //- Declare friendship with the const_iterator
        friend class const_iterator;


    // Constructors

        //- Construct given initial table size
        HashTableSH(const label size = 128);

        //- Construct from Istream
        HashTableSH(Istream&, const label size = 128);

        //- Construct as copy
        HashTableSH(const HashTableSH&);




    //- Destructor
    ~HashTableSH();


    // Member Functions

        // Access

            //- The size of the underlying table
            inline label capacity() const;

            //- Return number of elements in table
            inline label size() const;

            //- Return true if the hash table is empty
            inline bool empty() const;

            //- Return true if hashedEntry is found in table
            bool found(const long int&) const;

            //- Find and return an iterator set at the hashedEntry
            //  If not found iterator = end()
            iterator find(const long int&);

            //- Find and return an const_iterator set at the hashedEntry
            //  If not found iterator = end()
            const_iterator find(const long int&) const;

            //- Return the table of contents
           // List<long int> toc() const;

            //- Return the table of contents as a sorted list
            //List<long int> sortedToc() const;

            //- Print information
            Ostream& printInfo(Ostream&) const;


        // Edit

            //- Insert a new hashedEntry
            inline bool insert(const long int&, const int& newElmt);

            //- Assign a new hashedEntry, overwriting existing entries
            inline bool set(const long int&, const int& newElmt);
			
					bool set(const long int&, const int&,  bool , int &);

            //- Erase a hashedEntry specified by given iterator
            //  This invalidates the iterator until the next operator++
            bool erase(const iterator&);

            //- Erase a hashedEntry specified by the given key
            bool erase(const long int&);

            //- Remove entries given by the listed keys from this HashTableSH
            //  Return the number of elements removed
           // label erase(const UList<long int>&);

            //- Remove entries given by the given keys from this HashTableSH
            //  Return the number of elements removed.
            //  The parameter HashTableSH needs the same type of key, but the
            //  type of values held and the hashing function are arbitrary.
           // template<class AnyType, class AnyHash>
            label erase(const HashTableSH&);

            //- Resize the hash table for efficiency
            void resize(const label newSize);

            //- Clear all entries from table
            void clear();

            //- Clear the table entries and the table itself.
            //  Equivalent to clear() followed by resize(0)
            void clearStorage();

            //- Shrink the allocated table to approx. twice number of elements
            void shrink();

            //- Transfer the contents of the argument table into this table
            //  and annul the argument table.
            void transfer(HashTableSH&);



    // Member Operators

        //- Find and return a hashedEntry
        inline int& operator[](const long int&);

        //- Find and return a hashedEntry
        inline const int& operator[](const long int&) const;

        //- Find and return a hashedEntry, create it null if not present
        inline int& operator()(const long int&);

        //- Assignment
        void operator=(const HashTableSH&);

        //- Equality. Hash tables are equal if the keys and values are equal.
        //  Independent of table storage size and table order.
        bool operator==(const HashTableSH&) const;

        //- The opposite of the equality operation. Takes linear time.
        bool operator!=(const HashTableSH&) const;



    // STL type definitions

        //- Type of values the HashTableSH contains.
        typedef int value_type;

        //- Type that can be used for storing into HashTableSH::value_type
        //  objects.  This type is usually List::value_type&.
        typedef int& reference;

        //- Type that can be used for storing into constant
        //  HashTableSH::value_type objects.  This type is usually const
        //  HashTableSH::value_type&.
        typedef const int& const_reference;

        //- The type that can represent the size of a HashTableSH.
        typedef label size_type;


    // Iterators and helpers

        //- The iterator base for HashTableSH
        //  Note: data and functions are protected, to allow reuse by iterator
        //  and prevent most external usage.
        //  iterator and const_iterator have the same size, allowing
        //  us to reinterpret_cast between them (if desired)
        class iteratorBase
        {
            // Private Data

                //- Pointer to the HashTableSH for which this is an iterator
                //  This also lets us use the default bitwise copy/assignment
                HashTableSH* hashTable_;

                //- Current element
                hashedEntry* entryPtr_;

                //- Current hash index
                label hashIndex_;


        protected:

            // Constructors

                //- Construct null - equivalent to an 'end' position
                inline iteratorBase();

                //- Construct from hash table, moving to its 'begin' position
                inline explicit iteratorBase
                (
                    const HashTableSH* curHashTableSH
                );

                //- Construct from hash table, element and hash index
                inline iteratorBase
                (
                    const HashTableSH* curHashTableSH,
                    const hashedEntry* elmt,
                    const label hashIndex
                );


            // Protected Member Functions

                //- Increment to the next position
                inline void increment();

                //- Erase the HashTableSH element at the current position
                bool erase();

                //- Return non-const access to referenced object
                inline int& object();

                //- Return const access to referenced object
                inline const int& cobject() const;


        public:

            // Member operators

                // Access

                //- Return the long int corresponding to the iterator
                inline const long int& key() const;

                //- Compare hashedEntry element pointers
                inline bool operator==(const iteratorBase&) const;
                inline bool operator!=(const iteratorBase&) const;

                //- Compare hashedEntry to iteratorEnd pointers
                inline bool operator==(const iteratorEnd& unused) const;
                inline bool operator!=(const iteratorEnd& unused) const;
        };


        //- An STL-conforming iterator
        class iterator
        :
            public iteratorBase
        {
            friend class HashTableSH;

            // Private Member Functions

                //- Construct from hash table, moving to its 'begin' position
                inline explicit iterator
                (
                    HashTableSH* curHashTableSH
                );

                //- Construct from hash table, element and hash index
                inline iterator
                (
                    HashTableSH* curHashTableSH,
                    hashedEntry* elmt,
                    const label hashIndex
                );


        public:

            // Constructors

                //- Construct null (end iterator)
                inline iterator();

                //- Construct end iterator
                inline iterator(const iteratorEnd& unused);


            // Member operators

                //- Return referenced hash value
                inline int& operator*();
                inline int& operator()();

                //- Return referenced hash value
                inline const int& operator*() const;
                inline const int& operator()() const;

                inline iterator& operator++();
                inline iterator operator++(int);
        };

        //- Iterator set to the beginning of the HashTableSH
        inline iterator begin();


    // STL const_iterator

        //- An STL-conforming const_iterator
        class const_iterator
        :
            public iteratorBase
        {
            friend class HashTableSH;

            // Private Member Functions

                //- Construct from hash table, moving to its 'begin' position
                inline explicit const_iterator
                (
                    const HashTableSH* curHashTableSH
                );

                //- Construct from hash table, element and hash index
                inline const_iterator
                (
                    const HashTableSH* curHashTableSH,
                    const hashedEntry* elmt,
                    const label hashIndex
                );


        public:

            // Constructors

                //- Construct null (end iterator)
                inline const_iterator();

                //- Construct from iterator
                inline const_iterator(const iterator&);

                //- Construct end iterator
                inline const_iterator(const iteratorEnd& unused);


            // Member operators

                //- Return referenced hash value
                inline const int& operator*() const;
                inline const int& operator()() const;

                inline const_iterator& operator++();
                inline const_iterator operator++(int);

        };


        //- const_iterator set to the beginning of the HashTableSH
        inline const_iterator cbegin() const;

        //- const_iterator set to the beginning of the HashTableSH
        inline const_iterator begin() const;


    // IOstream Operator

        friend Istream& operator>> 
        (
            Istream&,
            HashTableSH&
        );

        friend Ostream& operator<< 
        (
            Ostream&,
            const HashTableSH&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    #include "HashTableSHI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//#ifndef NoHashTableSHC
//#define NoHashTableSHC
//#ifdef NoRepository
  // #include "HashTableSH.C"
//#endif
//#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
