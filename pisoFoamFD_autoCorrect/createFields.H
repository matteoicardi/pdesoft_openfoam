Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field lambda\n" << endl;
volScalarField lambda
(
    IOobject
    (
        "lambda",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero",dimless,0) // initialising to 0
);

Info<< "Reading field f\n" << endl;
volVectorField f
(
    IOobject
    (
        "f",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


#include "createPhi.H"


label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, piso.dict(), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());


singlePhaseTransportModel laminarTransport(U, phi);

autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, laminarTransport)
);

// * * * * *  Reading triSurface stl description of the geometry * * * * * * //        
IOdictionary fdDict
(
        IOobject
        (
            "fdDict",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
);

word fdName_( fdDict.lookup("fdName") );

dimensionedVector solidvel("vel", dimVelocity, fdDict.lookup("solidVel"));

bool use_stl(fdDict.lookup("useSTL"));

scalar alpha=readScalar(fdDict.lookup("alpha"));

scalar tolerance=readScalar(fdDict.lookup("tolerance"));
scalar tolerance_=max(1e-20,min(1e-2,tolerance));

if (use_stl)
{
triSurfaceMesh fdSurf_
(
    IOobject
    (
        fdName_+".stl",
        runTime.constant(),
        "triSurface",
        runTime,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

Info<< "Read Immersed Boundary triSurface" << endl;
fdSurf_.writeStats(Info);

triSurface fdTemp( fdSurf_ );
triSurfaceSearch fdTriSurfSearch( fdTemp );
boolList centresInside = fdTriSurfSearch.calcInside( mesh.cellCentres() );

forAll(mesh.C(),cellI)
{
    if ( centresInside[cellI] )
    {
        lambda[cellI] = 1.;
    }
}
}
