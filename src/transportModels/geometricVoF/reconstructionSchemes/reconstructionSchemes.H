/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019-2019 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2019-2019 DLR
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
    Foam::reconstructionSchemes

Description

    Original code supplied by Henning Scheufler, DLR (2019)

SourceFiles
    reconstructionSchemes.C

\*---------------------------------------------------------------------------*/

#ifndef reconstructionSchemes_H
#define reconstructionSchemes_H

#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"
#include "dimensionedScalar.H"
#include "autoPtr.H"
#include "IOdictionary.H"
#include "DynamicField.H"



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class reconstructionSchemes Declaration
\*---------------------------------------------------------------------------*/

class reconstructionSchemes
:
    public IOdictionary // will write the pointer in registry
    // inheritance from MeshObject(updateable) caused crashes with the polyhedralMesh
    // probably caused by a portation bug
{

private:

    // Private data

	    //- Coefficients dictionary
        dictionary reconstructionSchemesCoeffs_;


protected:

    // Protected data

        volScalarField& alpha1_;

        const surfaceScalarField& phi_;

        const volVectorField& U_;

        volVectorField normal_;

        volVectorField centre_;

        boolList interfaceCell_;

        DynamicField<label> interfaceLabels_;

        bool writeVTK_;

        Pair<label> timeIndexAndIter_;

        //- Write face points to .vtk files e.g. for visualisation in Paraview
        void isoFacesToFile
        (
            const DynamicList< List<point> >& isoFacePts,
            const word fileName,
            const word fileDir
        ) const;

        //- Write face points to .vtk files e.g. for visualisation in Paraview
        void isoFacesToFile
        (
            const List<word>& fieldName,
            const List<DynamicList<double> > &values,
            const DynamicList< List<point> >& faces,
            const word filNam,
            const word filDir
        ) const;

        bool alreadyReconstructed();


    // Private Member Functions

        //- No copy construct
        reconstructionSchemes(const reconstructionSchemes&) = delete;

        //- No bitwise assignment
        void operator=(const reconstructionSchemes&) = delete;


public:

	//- Return access to the model dictionary
	dictionary& modelDict();

	//- Return const access to the model dictionary
	const dictionary& modelDict() const;


public:

    //- Runtime type information
    TypeName("reconstructionSchemes");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            reconstructionSchemes,
            components,
            (
                volScalarField& alpha1,
                const surfaceScalarField& phi,
                const volVectorField& U,
                dictionary& dict
            ),
            (alpha1, phi,U,dict)
        );


    // Selectors

        //- Return a reference to the selected phaseChange model
        static autoPtr<reconstructionSchemes> New
        (
            volScalarField& alpha1,
            const surfaceScalarField& phi,
            const volVectorField& U,
            dictionary& dict
        );


    // Constructors

        //- Construct from components
        reconstructionSchemes
        (
            const word& type,
            volScalarField& alpha1,
            const surfaceScalarField& phi,
            const volVectorField& U,
            dictionary& dict
        );


    //- Destructor
    virtual ~reconstructionSchemes() = default;


    // Member Functions

        virtual void reconstruct() = 0;

        const volVectorField& normal() const
        {
            return normal_;
        }

        const volVectorField& centre() const
        {
            return centre_;
        }

        const boolList& interfaceCell() const
        {
            return interfaceCell_;
        }

        const DynamicField<label>& interfaceLabels() const
        {
            return interfaceLabels_;
        }

        virtual void mapAlphaField() const
        {
            // do nothing;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
