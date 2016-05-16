/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2016 OpenFOAM Foundation
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
    Foam::functionObjects::fieldValues::fieldValueDelta

Group
    grpFieldFunctionObjects

Description
    This function object provides applies an operation to the output of two
    fieldValue function objects.

    The operation is applied to all results of each fieldValue object.
    Accordingly, each object must generate the same number and type of results.

    Example of function object specification:
    \verbatim
    fieldValueDelta1
    {
        type            fieldValueDelta;
        libs ("libfieldFunctionObjects.so");
        operation       subtract;

        source1
        {
            ...
        }
        source2
        {
            ...
        }
    }
    \endverbatim

    \heading Function object usage
    \table
        Property     | Description             | Required    | Default value
        type         | type name: fieldValueDelta   | yes    |
    \endtable

    \linebreak
    The \c operation is one of:
    \plaintable
       add           | add
       subtract      | subtract
       min           | minimum
       max           | maximum
       average       | average
    \endplaintable

SeeAlso
    Foam::fieldValue

SourceFiles
    fieldValueDelta.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_fieldValueDelta_H
#define functionObjects_fieldValueDelta_H

#include "writeFiles.H"
#include "fieldValue.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
namespace fieldValues
{

/*---------------------------------------------------------------------------*\
                       Class fieldValueDelta Declaration
\*---------------------------------------------------------------------------*/

class fieldValueDelta
:
    public writeFiles
{
public:
        //- Operation type enumeration
        enum operationType
        {
            opAdd,
            opSubtract,
            opMin,
            opMax,
            opAverage
        };

        //- Operation type names
        static const NamedEnum<operationType, 5> operationTypeNames_;


private:

    // Private data

        //- Operation to apply to values
        operationType operation_;

        //- Field value source object 1
        autoPtr<fieldValue> source1Ptr_;

        //- Field value source object 2
        autoPtr<fieldValue> source2Ptr_;


    // Private Member Functions

        //- Templated function to apply the operation
        template<class Type>
        void apply
        (
            const word& resultType,
            const word& name1,
            const word& name2,
            const word& entryName1,
            const word& entryName2,
            bool& found
        );


protected:

    // Protected Member Functions

        //- Output file header information
        virtual void writeFileHeader(Ostream& os) const;


public:

    //- Run-time type information
    TypeName("fieldValueDelta");


    // Constructors

        //- Construct from Time and dictionary
        fieldValueDelta
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~fieldValueDelta();


    // Public Member Functions

        //- Read from dictionary
        virtual bool read(const dictionary&);

        //- Do nothing
        virtual bool execute(const bool postProcess = false);

        //- Calculate and write
        virtual bool write(const bool postProcess = false);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fieldValues
} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fieldValueDeltaTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
