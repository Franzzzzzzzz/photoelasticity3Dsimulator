



## Install
- Simplest is to start by installing [Miniconda](https://docs.anaconda.com/miniconda/install/)
- Install [Fenicsx](https://fenicsproject.org/download/)
- Install gmsh (in the appropriate conda environment) `conda install -c conda-forge gmsh python-gmsh`
- Install SDL2: `sudo apt install libsdl2-dev`
- Compile & build: `mkdir build && cd build`
- `cmake ..`

## Parameters
The configuration can be read from a json5 file directly provided on the commandline (`./photoelasticity config.json5`). The following key-values are recognised:
- `photoelastic constant`: 
  - Single number: constant for the monochrome photoelasticity
  - array of 3 numbers: constants for the reg, green, and blue wavelengths. 
  - arrya of 4 numbers: constants for the monochrome, red, green, and blue wavelengths, in that order. 
- `absorption`: double, default:0. absorbe the ray as well. Can be used to simulate X-ray for example. NOT IMPLEMENTED. 
- `color`: boolean, default: `false`. Indicates if the photoelasticity should be processed in RGB or with a single photoelastic value. 
- `pre-polarisation`: polarisation of the incident ray.
- `post-polarisation`: polariser for the ray just before hitting the detector. 

The polarisation can be given as a single number, corresponding to the angle from horizontal of the fast axis of a linear polariser, or one of the words "horizontal", "vertical", "right", "left", the latter two corresponding to circular polarisation. 

- `mesh`: filename of the mesh to use for the FEM. 

- `image size: [width, height]`: size of the image to produce, in pixels. 
- `pixel size`: physical size of the pixels of the image (in distance unit). 
- `display size: [width, height]`: size of the display, in pixels. The displayed imahe is interpolated from the ray calculation on the image size.  

- `sampling`: integer, default:100. Sampling rate for the ray propagation. 
- `strategy`: modify the internal calculation performed by the software. The following are currently defined: 
  - LINEAR_NEARESTNEIGHBOUR: equidistant sampling along the ray within the grain, each point is allocated the stress in the "nearest" tetrahedron, calculated from the center of mass of the tetrahedron.
  - LINEAR_TETRAHEDRON_INVERSION: equidistant sampling along the ray within the grain, however the intersection with each tetrahedra within the ray are actually calculated, and the stress used at each point is the one from the tetrahedron it belongs to.
  - LINEAR_TETRAHEDRON_MOLLERTRUMBORE: same as previous, but the ray-tetrahedron intersection is calculated using the Möller Trumbore algorithm, rather than the matrix inversion (my own algo). 
  - TETRAHEDRON_EXPONENTIAL_INVERSION: the intersection with tetrahedra is calculated, and the polarisation transformation matrix is calculated from the actual matrix exponentiation, rather than repeated exponentiation as with the previous methods.
  - TETRAHEDRON_EXPONENTIAL_MOLLERTRUMBORE: same as previous, but the ray-tetrahedron intersection is calculated using the Möller Trumbore algorithm, rather than the matrix inversion (my own algo). 

## Interactive mode shortcuts
Several key shortcuts are defined in interactive mode to modify on the fly the displayed image:
