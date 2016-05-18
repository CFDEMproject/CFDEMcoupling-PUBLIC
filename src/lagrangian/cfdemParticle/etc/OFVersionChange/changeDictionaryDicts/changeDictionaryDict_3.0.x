/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  3.0.x                                 |
|   \\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      changeDictionaryDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dictionaryReplacement
{

    turbulenceProperties 
    {
        simulationType  laminar;

        RAS
        {
            RASModel        laminar;

            turbulence      off;

            printCoeffs     on;
        }
    }
    couplingProperties
    {
        turbulenceModelType turbulenceProperties;
    }
}
// ************************************************************************* //
