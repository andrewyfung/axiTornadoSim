/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2020 ENERCON GmbH
    Copyright (C) 2018-2020 OpenCFD Ltd
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

#include "actuationDiskSource.H"
#include "fvMesh.H"
#include "fvMatrix.H"
#include "volFields.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class AlphaFieldType, class RhoFieldType>
void Foam::fv::actuationDiskSource::calc
(
    const AlphaFieldType& alpha,
    const RhoFieldType& rho,
    fvMatrix<vector>& eqn
)
{
    switch (forceMethod_)
    {
        case forceMethodType::FROUDE:
        {
            calcFroudeMethod(alpha, rho, eqn);
            break;
        }

        case forceMethodType::VARIABLE_SCALING:
        {
            calcVariableScalingMethod(alpha, rho, eqn);
            break;
        }

        case forceMethodType::TS_MAIN:
        {
            calcTSMain(alpha, rho, eqn);
            break;
        }

        case forceMethodType::TS_TURN:
        {
            calcTSTurn(alpha, rho, eqn);
            break;
        }

        case forceMethodType::POROUS:
        {
            calcPorous(alpha, rho, eqn);
            break;
        }

        default:
            break;
    }
}


template<class AlphaFieldType, class RhoFieldType>
void Foam::fv::actuationDiskSource::calcFroudeMethod
(
    const AlphaFieldType& alpha,
    const RhoFieldType& rho,
    fvMatrix<vector>& eqn
)
{
    const vectorField& U = eqn.psi();
    vectorField& Usource = eqn.source();
    const scalarField& cellsV = mesh_.V();

    // Compute upstream U and rho, spatial-averaged over monitor-region
    vector Uref(Zero);
    scalar rhoRef = 0.0;
    label szMonitorCells = monitorCells_.size();

    for (const auto& celli : monitorCells_)
    {
        Uref += U[celli];
        rhoRef = rhoRef + rho[celli];
    }
    reduce(Uref, sumOp<vector>());
    reduce(rhoRef, sumOp<scalar>());
    reduce(szMonitorCells, sumOp<label>());

    if (szMonitorCells == 0)
    {
        FatalErrorInFunction
            << "No cell is available for incoming velocity monitoring."
            << exit(FatalError);
    }

    Uref /= szMonitorCells;
    rhoRef /= szMonitorCells;

    const scalar Ct = sink_*UvsCtPtr_->value(mag(Uref));
    const scalar Cp = sink_*UvsCpPtr_->value(mag(Uref));

    if (Cp <= VSMALL || Ct <= VSMALL)
    {
        FatalErrorInFunction
           << "Cp and Ct must be greater than zero." << nl
           << "Cp = " << Cp << ", Ct = " << Ct
           << exit(FatalIOError);
    }

    // (BJSB:Eq. 3.9)
    const scalar a = 1.0 - Cp/Ct;
    const scalar T = 2.0*rhoRef*diskArea_*magSqr(Uref & diskDir_)*a*(1 - a);

    for (const label celli : cells_)
    {
        Usource[celli] += ((cellsV[celli]/V())*T)*diskDir_;
    }

    if
    (
        mesh_.time().timeOutputValue() >= writeFileStart_
     && mesh_.time().timeOutputValue() <= writeFileEnd_
    )
    {
        Ostream& os = file();
        writeCurrentTime(os);

        os  << Uref << tab << Cp << tab << Ct << tab << a << tab << T
            << endl;
    }
}


template<class AlphaFieldType, class RhoFieldType>
void Foam::fv::actuationDiskSource::calcVariableScalingMethod
(
    const AlphaFieldType& alpha,
    const RhoFieldType& rho,
    fvMatrix<vector>& eqn
)
{
    const vectorField& U = eqn.psi();
    vectorField& Usource = eqn.source();
    const scalarField& cellsV = mesh_.V();

    // Monitor and average monitor-region U and rho
    vector Uref(Zero);
    scalar rhoRef = 0.0;
    label szMonitorCells = monitorCells_.size();

    for (const auto& celli : monitorCells_)
    {
        Uref += U[celli];
        rhoRef = rhoRef + rho[celli];
    }
    reduce(Uref, sumOp<vector>());
    reduce(rhoRef, sumOp<scalar>());
    reduce(szMonitorCells, sumOp<label>());

    if (szMonitorCells == 0)
    {
        FatalErrorInFunction
            << "No cell is available for incoming velocity monitoring."
            << exit(FatalError);
    }

    Uref /= szMonitorCells;
    const scalar magUref = mag(Uref);
    rhoRef /= szMonitorCells;

    // Monitor and average U and rho on actuator disk
    vector Udisk(Zero);
    scalar rhoDisk = 0.0;
    scalar totalV = 0.0;

    for (const auto& celli : cells_)
    {
        Udisk += U[celli]*cellsV[celli];
        rhoDisk += rho[celli]*cellsV[celli];
        totalV += cellsV[celli];
    }
    reduce(Udisk, sumOp<vector>());
    reduce(rhoDisk, sumOp<scalar>());
    reduce(totalV, sumOp<scalar>());

    if (totalV < SMALL)
    {
        FatalErrorInFunction
            << "No cell in the actuator disk."
            << exit(FatalError);
    }

    Udisk /= totalV;
    const scalar magUdisk = mag(Udisk);
    rhoDisk /= totalV;

    if (mag(Udisk) < SMALL)
    {
        FatalErrorInFunction
            << "Velocity spatial-averaged on actuator disk is zero." << nl
            << "Please check if the initial U field is zero."
            << exit(FatalError);
    }

    // Interpolated thrust/power coeffs from power/thrust curves
    const scalar Ct = sink_*UvsCtPtr_->value(magUref);
    const scalar Cp = sink_*UvsCpPtr_->value(magUref);

    if (Cp <= VSMALL || Ct <= VSMALL)
    {
        FatalErrorInFunction
           << "Cp and Ct must be greater than zero." << nl
           << "Cp = " << Cp << ", Ct = " << Ct
           << exit(FatalIOError);
    }

    // Calibrated thrust/power coeffs from power/thrust curves (LSRMTK:Eq. 6)
    const scalar CtStar = Ct*sqr(magUref/magUdisk);
    const scalar CpStar = Cp*pow3(magUref/magUdisk);

    // Compute calibrated thrust/power (LSRMTK:Eq. 5)
    const scalar T = 0.5*rhoRef*diskArea_*magSqr(Udisk & diskDir_)*CtStar;
    const scalar P = 0.5*rhoRef*diskArea_*pow3(mag(Udisk & diskDir_))*CpStar;

    for (const label celli : cells_)
    {
        Usource[celli] += (cellsV[celli]/totalV*T)*diskDir_;
    }

    if
    (
        mesh_.time().timeOutputValue() >= writeFileStart_
     && mesh_.time().timeOutputValue() <= writeFileEnd_
    )
    {
        Ostream& os = file();
        writeCurrentTime(os);

        os  << Uref << tab << Cp << tab << Ct << tab
            << Udisk << tab << CpStar << tab << CtStar << tab << T << tab << P
            << endl;
    }
}



template<class AlphaFieldType, class RhoFieldType>
void Foam::fv::actuationDiskSource::calcTSMain
(
    const AlphaFieldType& alpha,
    const RhoFieldType& rho,
    fvMatrix<vector>& eqn
)
{
    const vectorField& U = eqn.psi();
    vectorField& Usource = eqn.source();
    const scalarField& cellsV = mesh_.V();
    //const volVectorField& cellsC = mesh_.C();

    // Monitor and average U and rho on actuator disk
    vector Udisk(Zero);
    scalar rhoDisk = 0.0;
    scalar totalV = 0.0;
    for (const auto& celli : cells_)
    {
        Udisk += U[celli]*cellsV[celli];
        rhoDisk += rho[celli]*cellsV[celli];
        totalV += cellsV[celli];
    }
    reduce(Udisk, sumOp<vector>());
    reduce(rhoDisk, sumOp<scalar>());
    reduce(totalV, sumOp<scalar>());

    if (totalV < SMALL)
    {
        FatalErrorInFunction
            << "No cell in the actuator disk."
            << exit(FatalError);
    }

    Udisk /= totalV;
    rhoDisk /= totalV;

    // Compute calibrated thrust/power (LSRMTK:Eq. 5)
    const scalar T = diskArea_;

    Ostream& os = file();

    for (const label celli : cells_)
    {
        // Read into blade-element method before trying to use this again.
        //const scalar fanRadius = sqr(cellsC[celli].x()) + sqr(cellsC[celli].z());
        //const scalar tanFanBeta = 0.01; //Tangent of fan blade angle
        //vector dir(-cellsC[celli].z()*tanFanBeta*(fanRadius/0.87), -1*(1-tanFanBeta), cellsC[celli].x()*tanFanBeta*(fanRadius/0.87));

        vector dir(0, -1, 0);

        Usource[celli] += T*(cellsV[celli]/totalV)*dir;
    }

    if
    (
        mesh_.time().timeOutputValue() >= writeFileStart_
     && mesh_.time().timeOutputValue() <= writeFileEnd_
    )
    {
        writeCurrentTime(os);

        os  << Udisk << T << tab << endl;
    }
}



template<class AlphaFieldType, class RhoFieldType>
void Foam::fv::actuationDiskSource::calcTSTurn
(
    const AlphaFieldType& alpha,
    const RhoFieldType& rho,
    fvMatrix<vector>& eqn
)
{
    const vectorField& U = eqn.psi();
    vectorField& Usource = eqn.source();
    const scalarField& cellsV = mesh_.V();
    const volVectorField& cellsC = mesh_.C();


    Ostream& os = file();

    // Monitor and average monitor-region U and rho
    vector Uref(Zero); //---
    scalar flowAngle = 0.0;
    label szMonitorCells = monitorCells_.size();

    for (const auto& celli : monitorCells_)
    {
        Uref += U[celli]; //---

        // scalar r = sqrt(pow(cellsC[celli].x(),2) + pow(cellsC[celli].z(),2));

        scalar angle = atan(cellsC[celli].z()/(cellsC[celli].x()+0.000001));
        
        scalar Urad = mag(U[celli].x()*cos(angle)) + mag(U[celli].z()*sin(angle));
        scalar Utan = mag(U[celli].z()*cos(angle)) - mag(U[celli].x()*sin(angle));

        flowAngle = flowAngle + mag(atan(Utan/(Urad+0.000001)));

        // os << cellsC[celli].x() << tab << cellsC[celli].y() << tab << cellsC[celli].z() << tab << angle << endl;
        // os << tab << U[celli].x() << tab << U[celli].z() << tab << Utan << tab << Uy << endl;

    }
    reduce(Uref, sumOp<vector>()); //---
    reduce(flowAngle, sumOp<scalar>());
    reduce(szMonitorCells, sumOp<label>());

    Uref /= szMonitorCells;
    flowAngle /= szMonitorCells;
    flowAngle *= 180/3.14159;

    // Monitor and average U and rho on actuator disk
    vector Udisk(Zero);
    scalar rhoDisk = 0.0;
    scalar totalV = 0.0;
    for (const auto& celli : cells_)
    {
        Udisk += U[celli]*cellsV[celli];
        rhoDisk += rho[celli]*cellsV[celli];
        totalV += cellsV[celli];
    }
    reduce(Udisk, sumOp<vector>());
    reduce(rhoDisk, sumOp<scalar>());
    reduce(totalV, sumOp<scalar>());

    if (totalV < SMALL)
    {
        FatalErrorInFunction
            << "No cell in the actuator disk."
            << exit(FatalError);
    }

    Udisk /= totalV;
    rhoDisk /= totalV;


    // Compute calibrated thrust/power
    if (diskDir_[1] > 0.5) {
        const scalar targetAngle = UvsCtPtr_->value(1);
        diskArea_ += UvsCpPtr_->value(1)*(targetAngle - flowAngle);
    }
    const scalar T = diskArea_;

    for (const label celli : cells_)
    {
        scalar r = sqrt(sqr(cellsC[celli].z()) + sqr(cellsC[celli].x()));
        vector rotDirection(-cellsC[celli].z()/r, 0, cellsC[celli].x()/r);

        Usource[celli] += T*(cellsV[celli]/totalV)*rotDirection;
    }

    if
    (
        mesh_.time().timeOutputValue() >= writeFileStart_
     && mesh_.time().timeOutputValue() <= writeFileEnd_
    )
    {
        writeCurrentTime(os);

        os << Uref << tab << flowAngle << tab << diskArea_ << tab << Udisk << tab << endl;
    }
}



template<class AlphaFieldType, class RhoFieldType>
void Foam::fv::actuationDiskSource::calcPorous
(
    const AlphaFieldType& alpha,
    const RhoFieldType& rho,
    fvMatrix<vector>& eqn
)
{
    const vectorField& U = eqn.psi();
    vectorField& Usource = eqn.source();
    const scalarField& cellsV = mesh_.V();
    //const volVectorField& cellsC = mesh_.C();

    // Monitor and average U and rho on actuator disk
    vector Udisk(Zero);
    scalar rhoDisk = 0.0;
    scalar totalV = 0.0;
    for (const auto& celli : cells_)
    {
        Udisk += U[celli]*cellsV[celli];
        rhoDisk += rho[celli]*cellsV[celli];
        totalV += cellsV[celli];
    }
    reduce(Udisk, sumOp<vector>());
    reduce(rhoDisk, sumOp<scalar>());
    reduce(totalV, sumOp<scalar>());

    if (totalV < SMALL)
    {
        FatalErrorInFunction
            << "No cell in the actuator disk."
            << exit(FatalError);
    }

    Udisk /= totalV;

    Ostream& os = file();

    const scalar C2_axial = 0.1;
    const scalar C2_tan = 1000;
    for (const label celli : cells_)
    {
        vector Ucell = U[celli];
        scalar Umag = mag(Ucell);
        scalar Ux = Ucell.x();
        scalar Uy = Ucell.y();
        scalar Uz = Ucell.z();
        //scalar r = sqrt(sqr(cellsC[celli].z()) + sqr(cellsC[celli].x()));
        //vector rotDirection(-cellsC[celli].z()/r, 0, cellsC[celli].x()/r);

        scalar Tmag = 1;

        vector Tvec(C2_tan*Umag*Ux, C2_axial*Umag*Uy, C2_tan*Umag*Uz);

        Usource[celli] += Tmag*Tvec*cellsV[celli];
        Info << "test" << endl;
        // Unsure if volume div is needed?
        //Usource[celli] += T*(cellsV[celli]/totalV);
    }

    if
    (
        mesh_.time().timeOutputValue() >= writeFileStart_
     && mesh_.time().timeOutputValue() <= writeFileEnd_
    )
    {
        writeCurrentTime(os);

        os  << Udisk << endl;
    }
}



// ************************************************************************* //
