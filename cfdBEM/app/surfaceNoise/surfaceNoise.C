/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

Application
    surfaceNoise

Description
    Utility to perform noise analysis of pressure data in points of surface using the libFFT
    library.

    Control settings are read from the $FOAM_CASE/system/surfaceNoiseDict dictionary,
    or user-specified dictionary using the -dict option. Pressure data prepared by
    surfacePressureSampling utility and read using a simple file reader:

    \heading Usage

    \verbatim

    outputFormat    gmsh; //optionally

    inputFileName   "fileName.txt";

    maxFrequency    1000;

    \endverbatim


    Output:
    - FFT of the pressure data

SeeAlso
    FoamFftwDriver.H

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "timeSelector.H"
#include "Time.H"
#include "Pair.H"

#include "FoamFftwDriver.H"
#include "FileInterface/FileInterface.H"
#include "Interpolator/Interpolator.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //



// just return Re and Im parts of complex pressure amplitude; frequency is out of function
Foam::autoPtr<Foam::List<Foam::complex> > fft(const List<scalar> &p, scalar lastTime)
{
    autoPtr<List<complex> > fft_resPtr;

    if ( (p.size() > 0) )
    {
        FoamFftwDriver fftw (p, lastTime);
        
        fft_resPtr = fftw.simpleComplexForwardTransform();
    }
    
    return fft_resPtr;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions(true, true);
    argList::noParallel();

    #include "addDictOption.H"
    #include "setRootCase.H"
    #include "createTime.H"

    #include "createFields.H"

    Info<< "Reading data file..." << endl;

    fileName pFileName = "pressureData/" + inputFileName;

    FileInterface reader(pFileName);

    autoPtr<List<List<scalar> > > ptData = reader.read();

    //------------------------
    //prepare place for sampling

    //number of sampled (init) times
    label sampledTimePoints = ptData().size();

    List<scalar> tSampled(sampledTimePoints,0.0);
    //tSampled.resize(sampledTimePoints);

    forAll(tSampled,tI)
    {
        tSampled[tI] = ptData()[tI][0];
    }

    scalar lastTime = tSampled[sampledTimePoints-1];

    // pressure data
    List<scalar> pSampled(sampledTimePoints,0.0); // just buffer for pressure data in one point
    //pSampled.resize(sampledTimePoints);

    //------------------------
    //calculate the interpolation time interval

    Info << "Prepare to fft..." << endl;

    label nFftTimePoints = label(2*maxFreq*lastTime); //2*... because spectrum is symmetric
    scalar deltaIntT = lastTime / (scalar(nFftTimePoints) );
    label freqPoints = ceil(0.5*scalar(nFftTimePoints)) + 1;

    //pInt.resize(nFftTimePoints);


    //------------------------
    //fft for all points on surface point-by-point

    List<scalar> freq(freqPoints,0.0);

    forAll(freq,i)
    {
        freq[i] = scalar(i)/lastTime;
    }

    List<List<complex > > allComplexPressures(freqPoints); // buffer for all complex pressures
    //allComplexPressures.resize(freqPoints);   

    autoPtr<List<complex > > fftOnePoint;

    label nSurfPoints = ptData()[0].size();

    Info << "Calculate interpolation and FFT..." << endl;

    for (label i = 1; i <= nSurfPoints; ++i) //loop starts from 1 because time value is placed on 0 position
    {                    
        forAll(pSampled,pI)
        {
            pSampled[pI] = ptData()[pI][i];
        }

        Interpolator interpolator;

        List<scalar> pInt = interpolator.interpolate(tSampled,pSampled,nFftTimePoints,deltaIntT);

        // calculate fft
        fftOnePoint = fft(pInt, lastTime);

        //write to buffer
        forAll(allComplexPressures,freqI)
        {
            allComplexPressures[freqI].resize(i);
            allComplexPressures[freqI][i-1] = fftOnePoint()[freqI];
        }
    }

    //------------------------

    Info<< "Write gmsh files..." << endl;

    fileName outputDirectory = inputFileName.lessExt();    

    fileName controlFile = "complexPressureData/" + outputDirectory + "/" + "Spectrum.dat";

    autoPtr<OFstream> controlFilePtr;

    if (controlFilePtr.empty())
    {
        // Open new file at start up
        controlFilePtr.reset
        (
            new OFstream
            (
                controlFile
            )
        );
    }

    forAll(freq, freqI)
    {
        word freqValue = name(freq[freqI]);

        fileName freqDirectory = "complexPressureData/" + outputDirectory + "/" + freqValue;

        //Info << freqDirectory << endl;

        mkDir(freqDirectory);

        //write appropriate frequency
        //fileName freqFile  = freqDirectory + "/frequency";

        //autoPtr<OFstream> freqFilePtr;

        //if (freqFilePtr.empty())
        //{
        // Open new file at start up
        //    freqFilePtr.reset
        //    (
        //        new OFstream
        //        (
        //            freqFile
        //        )
        //    );
        //}

        //freqFilePtr() << freq[freqI];

        //write data to msh file
        if (outputFormat == "gmsh")
        {
            fileName dataFile = freqDirectory + "/nodeData.msh";

            FileInterface outputFile(dataFile);

            outputFile.write(freq[freqI],allComplexPressures[freqI]);
        }

        Info << "---------------\n";
        Info << "Frequency: " << freq[freqI] << nl;

        scalar maxAbs = 0.0;

        forAll(allComplexPressures[freqI],k)
        {
            complex currentCA = allComplexPressures[freqI][k];
            scalar currentAbs = Foam::sqrt( pow(currentCA.Re(),2) + pow(currentCA.Im(),2) );

            if (maxAbs < currentAbs)
            {
                maxAbs = currentAbs;
            }
        }

        Info << "Max abs: " << maxAbs << nl;

        controlFilePtr() << freq[freqI] << ' ' << maxAbs << nl;
    }

    return 0;
}


// ************************************************************************* //
