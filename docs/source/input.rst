GraGLeS2D is a OpenMP-parallelized computer program for simulating the evolution of anisotropic grain growth phenomena in polycrystal that makes use of the level-set method. It is free to download from:: 

  https://github.com/cmiessen/GraGLeS2D

  https://github.com/GraGLeS/GraGLeS2D

GraGLeS 2D grain growth data format
===================================

For each time step at which a snapshot of the entire microstructure was taken, by default GraGLeS2D writes a plain ASCII file with the following syntax **Texture_<no>.txt**,
where **<no>** is positiv and refers to the simulation step ID at which the snapshot was taken. 

Thus, a collection of *N* snapshots contains all pieces of information to extract the temporal trajectory of individual grains.
Each such file, which is headerless and space-separated, contains always all still existent grains in the network, their properties and their neighbors.
Specifically, the file layout reads as follows:

* **Identification ID** - a positive grain label
* **Number of faces** - how many neighbors does the grain possess
* **BoundaryFlag** - is the grain located in contact with the simulation domain boundary (1 - yes, 0 - no)
* **Area** - size normalized to the size (1.0) of the unit square simulation domain
* **AreaChange** - size difference compared to previous time step 
* **Perimeter** - circumference
* **TotalBoundaryEnergy** - face length times energy for each face
* **OrientationBunge1** - Bunge convention Euler angle \phi_1
* **OrientationBunge2** - Bunge convention Euler angle \PHI
* **OrientationBunge3** - Bunge convention Euler angle \phi_2
* **ID FaceLength** - a list of exactly *Number of faces* data pairs indicating all neighbors with their ID and the length of the adjointing face

An example file is included in the *example/* folder.
