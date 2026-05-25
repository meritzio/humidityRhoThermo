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
    Foam::humidityRhoThermo

Description
    Base-class for fluid thermodynamic properties based on density.

See also
    Foam::fluidThermo

SourceFiles
    humidityRhoThermo.C

\*---------------------------------------------------------------------------*/

#ifndef humidityRhoThermo_H
#define humidityRhoThermo_H

#include "HumidityRhoThermo.H"
#include "pureThermo.H"
#include "fluidThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class humidityRhoThermo Declaration
\*---------------------------------------------------------------------------*/

class humidityRhoThermo
:
    virtual public fluidThermo
{
public:

    // Public Classes

        //- Forward declare the implementation class
        class implementation;

        //- Forward declare the composite class
        class composite;

    // Public Typedefs

        //- The derived type
        template<class MixtureType>
        using DerivedThermoType =
            HumidityRhoThermo<BasicThermo<MixtureType, composite>>;

        //- The derived name
        static word derivedThermoName()
        {
            return "heHumidityRhoThermo";
        }

    //- Runtime type information
    TypeName("humidityRhoThermo");

    //- Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        humidityRhoThermo,
        fvMesh,
        (const fvMesh& mesh, const word& phaseName),
        (mesh, phaseName)
    );

    // Selectors

        //- Standard selection based on fvMesh
        static autoPtr<humidityRhoThermo> New
        (
            const fvMesh& mesh,
            const word& phaseName=word::null
        );


    //- Destructor
    virtual ~humidityRhoThermo();


    // Member Functions

        // Fields derived from thermodynamic state variables

            //- Density [kg/m^3]
            virtual tmp<volScalarField> rho() const = 0;

            //- Density for patch [kg/m^3]
            virtual tmp<scalarField> rho(const label patchi) const = 0;

            //- Return the thermodynamic density field [kg/m^3]
            //  This is used by solvers which create a separate continuity rho
            virtual tmp<volScalarField> renameRho();

            //- Return non-const access to the local density field [kg/m^3]
            virtual volScalarField& rho() = 0;

            //- Add the given density correction to the density field.
            //  Used to update the density field following pressure solution.
            virtual void correctRho(const volScalarField& deltaRho);
};


/*---------------------------------------------------------------------------*\
                  Class humidityRhoThermo::implementation Declaration
\*---------------------------------------------------------------------------*/

class humidityRhoThermo::implementation
:
    virtual public humidityRhoThermo
{
protected:

    // Protected data

        //- Density field [kg/m^3]
        //  Named 'rho' to avoid (potential) conflict with solver density
        volScalarField rho_;
    
        //- Relative humidity [-]
        volScalarField relHum_;

        //- Water mass [kg]
        volScalarField waterMass_;

        //- Water vapor content [kg(water)/m^3(Air)]
        volScalarField waterVapor_;

        //- Max water vapor possible [kg(water)/m^3(Air)] for stabilization
        volScalarField maxWaterVapor_;

        //- Specific humidity [kg/kg] (water/air)
        volScalarField specificHumidity_;

        //- Max specific humidity [kg/kg] (water/air)
        volScalarField maxSpecificHumidity_;

        //- Saturation pressure of water [Pa]
        volScalarField pSatH2O_;

        //- Partial pressure of water [Pa]
        volScalarField partialPressureH2O_;

        //- Partial pressure calculation method
        mutable word method_;

        //- Initialize the specificHumidity using the relHum field
        bool initWithRelHumidity_;

        // Turbulence fields

        volScalarField muEff_;

public:

    // Constructors

        //- Construct from dictionary, mesh  and phase name
        implementation(const dictionary&, const fvMesh&, const word& phaseName);

        //- Disallow default bitwise copy construction
        implementation(const implementation&) = delete;


    //- Destructor
    virtual ~implementation();


    // Member Functions

        // Fields derived from thermodynamic state variables

            //- Density [kg/m^3]
            virtual tmp<volScalarField> rho() const;

            //- Density for patch [kg/m^3]
            virtual tmp<scalarField> rho(const label patchi) const;

            //- Return non-const access to the local density field [kg/m^3]
            virtual volScalarField& rho();

        // Read

            virtual void readMethod();

            virtual void readOrInitSpecificHumidity();
            
    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const implementation&) = delete;
};


/*---------------------------------------------------------------------------*\
                 Class humidityRhoThermo::composite Declaration
\*---------------------------------------------------------------------------*/



class humidityRhoThermo::composite
:
    public basicThermo::implementation,
    public pureThermo,
    public fluidThermo::implementation,
    public humidityRhoThermo::implementation
{
public:

    // Constructors

        //- Construct from dictionary, mesh and phase name
        template<class MixtureType>
        composite
        (
             const dictionary& dict,
             const MixtureType& mixture,
             const fvMesh& mesh,
             const word& phaseName
        )
        :
            basicThermo::implementation(dict, mesh, phaseName),
            fluidThermo::implementation(dict, mesh, phaseName),
            humidityRhoThermo::implementation(dict, mesh, phaseName)
        {}
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
