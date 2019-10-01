# Intergration VoFLib

* implicitFunction -> src/OpenFoam/primitives/functions
* VoFLibCore -> src/transportModels/geometricVoF
* interface -> src/sampling/surfaces
* markRegion -> src/transportModels/geometricVoF
* compressibleInterIsoFoam -> apps/solvers/multiphase/
* mod(interIsoFoam)
* mod(setAlphaField) -> should be compatible with the old version without changes of the dictionaries
* zoneDistribute -> src/finiteVolume/fvMesh
* leastSquare -> src/finiteVolume/fvMatrices/solvers
* RDF -> src/transportModels/geometricVoF
* added tests -> applications/test
    * setAlphaField
    * leastSquareGrad
    * multiDimPolyFitter
    * reconstructedDistanceFunction
    * zoneDistribute

## TODO

* update interIsoFoam tutorials reconstructionScheme and vof2IsoTol are missing
* setAlphaField createdSurface to VTK needs rework
* update Header

## Discussion

* renamed isoFaceTol to vof2IsoTol -> plic and isoAlpha in the same framework? deprecated keyword 

