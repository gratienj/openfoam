# Intergration VoFLib

* implicitFunction -> src/OpenFoam/primitives/functions
* VoFLibCore -> src/finiteVolume/geometricVoF
* interface -> src/sampling/surfaces
* markRegion -> functionObjects
* compressibleInterIsoFoam -> apps/solvers/multiphase/
* mod(interIsoFoam)
* mod(setAlphaField) -> should be compatible with the old version without changes of the dictionaries
* exchangeField -> src/finiteVolume/extendedStencil
* leastSquare -> src/finiteVolume/finiteVolumen
* RDF -> 


## TODO

* update interIsoFoam tutorials reconstructionScheme and vof2IsoTol are missing
* setAlphaField createdSurface to VTK needs rework
* update Header

## Discussion

* renamed isoFaceTol to vof2IsoTol -> plic and isoAlpha in the same framework? deprecated keyword 

