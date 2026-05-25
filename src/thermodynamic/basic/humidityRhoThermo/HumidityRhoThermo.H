/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2025 OpenFOAM Foundation
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
    Foam::HumidityRhoThermo

Description
    Energy for a mixture based on density

SourceFiles
    HumidityRhoThermo.C

\*---------------------------------------------------------------------------*/

#ifndef HumidityRhoThermo_H
#define HumidityRhoThermo_H

#include "BasicThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class HumidityRhoThermo Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo>
class HumidityRhoThermo
:
    public Thermo
{
    // Private Member Functions

        //- Calculate the thermo variables
        void calculate();


public:

    // Constructors

        //- Construct from mesh and phase name
        HumidityRhoThermo
        (
            const fvMesh&,
            const word& phaseName
        );

        //- Disallow default bitwise copy construction
        HumidityRhoThermo(const HumidityRhoThermo<Thermo>&) = delete;


    //- Destructor
    virtual ~HumidityRhoThermo();


    // Member Functions

        //- Initialize the thermo:specificHumidity using the thermo:relHum field
        virtual void initialize();

        //- Update properties
        virtual void correct();

        //- Calculate the saturation pressure of H2O [Pa]
        virtual void pSatH2O();

        //- Calculate partial pressure of water saturation pressure [Pa]
        //  based from the relative humidity (only for initialization)
        virtual void partialPressureH2OFromRelHumidity();

        //- Calculate partial pressure of water saturation pressure [Pa]
        virtual void partialPressureH2O();

        //- Calculate the relative humidity [-]
        virtual void relHumidity();

        //- Calculate the water vapor content [kg(water)/m^3(air)]
        virtual void waterVapor();

        //- Calculate the maximum specific humidity (for stabilization)
        virtual void maxSpecificHumidity();

        //- Calculate the water mass [kg(water)]
        virtual void waterMass();

        //- Conserved transport equation for the specific humidity [kg/kg]
        virtual void specificHumidityTransport();

        //- Initialize the specific humidity field by using the relative
        //  humidity field
        virtual void initialSpecificHumidityFromRelHumidity();

        //- Update density field
        virtual void updateRho(volScalarField&);

        //- Limit the specific humidity to the max range
        virtual void limitMax();
        
    // Member Operators

        //- Disallow default bitwise assignment
        void operator=
        (
            const HumidityRhoThermo<Thermo>&
        ) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "HumidityRhoThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
