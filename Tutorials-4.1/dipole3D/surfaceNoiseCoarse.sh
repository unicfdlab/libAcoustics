#!/bin/bash
#
# run surfaceNoise utility
#
cat system/surfaceNoiseDict_coarse > system/surfaceNoiseDict
surfaceNoise

surface="triSurfaceSamplingCoarse"
freqs=`ls complexPressureData/"$surface"_surfacePressure`

for f in $freqs
do
    echo "Processing frequency $f"
    cat surfaceGeometryData/$surface.msh > complexPressureData/"$surface"_surfacePressure/$f/dataWithMesh.msh
    cat complexPressureData/"$surface"_surfacePressure/$f/nodeData.msh >> complexPressureData/"$surface"_surfacePressure/$f/dataWithMesh.msh
done

#
#END-OF-FILE
#

