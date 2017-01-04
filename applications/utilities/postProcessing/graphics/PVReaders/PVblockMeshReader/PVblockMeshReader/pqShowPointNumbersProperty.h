/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenCFD Ltd.
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
    pqShowPointNumbersProperty

Description
    Custom UI handling of show-points (ParaView blockMesh reader)

SourceFiles
    pqShowPointNumbersProperty.cxx

\*---------------------------------------------------------------------------*/
#ifndef pqShowPointNumbersProperty_h
#define pqShowPointNumbersProperty_h

#include "pqPropertyWidget.h"

// Forward declarations (ParaView)
class vtkSMIntVectorProperty;


/*---------------------------------------------------------------------------*\
                 Class pqShowPointNumbersProperty Declaration
\*---------------------------------------------------------------------------*/

class pqShowPointNumbersProperty
:
    public pqPropertyWidget
{
    Q_OBJECT;
    typedef pqPropertyWidget Superclass;

    // Private data

        //- Show Point Numbers (bool property)
        vtkSMIntVectorProperty* showPointNumbers_;


protected slots:

    // Protected Member Functions

    //- Sync property with changed checkbox state, update rendered view(s)
    void showPointNumbers(bool checked);


public:

    //- Construct from components
    pqShowPointNumbersProperty
    (
        vtkSMProxy* proxy,
        vtkSMProperty* prop,
        QWidget* parent = nullptr
    );


    //- Destructor
    virtual ~pqShowPointNumbersProperty();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
