OpenAccessSPC
=============

Code for creating and using the SPC model of a plenoptic camera with known geometry.

The current implementation can be run with both MATLAB and GNU Octave.  It will 
be rewritten for Python, hopefully soon.


Generating the SPC-model
-------------------------------------------------------------------------------

The functionality for generating an SPC model for a camera system is 
implemented in the following files:

    * Generate_LCm.m
    * Generate_LCm_decoupled.m
    * Generate_Aperture.m
    * Aperture_m.m
    * Generate_Lens.m
    * Lens_m.m

Each file contains a function with the same name as the file.

The first function, `Generate_LCm.m` generates the initial set of light cones
for the image sensor pixels. The angular span of the generated light cones is 
defined with respect to the given light-acceptance-angle of the sensor pixel 
as a physical property of the image sensor pixel. More info about the 
parameters can be found in the header of the file. 

The function `Generate_LCm_decoupled()`, contains the functionality for 
generating the initial set of light cones on the sensor pixel, and assumes 
that all lenslets in a micro-lens array are decoupled from their neighbouring 
lenslets.  This means that each lenslet affects only the pixels behind that 
particular lenslet, see e.g. [SPIE13] for details.
More info about the parameters can be found in the header of the file.

The function `Generate_Aperture()` generates an aperture from given values.  
Then `Aperture_m()` applies that aperture to a set of light cones.  In this 
implementation the aperture describes the part of the lens that light can pass 
through, as otherwise the lens is assumed to be infinitely wide.
More info about the parameters can be found in the header of the file.

The function `Generate_Lens()` is then used to generate a lens or a lens array 
from given values.  Then the function `Lens_m()` is used to apply this lens or 
lens-array to the light cones generated above and trimmed by the aperture.
More info about the parameters can be found in the header of the file.

The use of the functions just described are illustrated in a set of example 
files:

    * PC2.m
    * Lytro.m
    * R29.m

All implementations are given as 1-dimensional optimizations.

In `PC2.m` I use the set of the previously defined functions to generate 
the SPC model of a plenoptic camera. An image sensor and a microlens-array 
with known geometry and spacing build the structure of the plenoptic camera.

In `Lytro.m` I've implemented the original Lytro camera as an example of  
plenoptic camera implementation. 
More info about the parameters can be found in the header of the file.
See [Lytro] for more details.

`R29.m` is my implementation of the Raytrix R29 plenoptic camera, which is used 
to generate the results in [Raytrix].


Extracting camera properties
-------------------------------------------------------------------------------

So far the special resolution has been the interesting property.  Here I look 
at the spatial resolution profile through depth, or, as I call it, the lateral 
resolution.

The following scripts are described:

    * Base_m.m
    * ResDistPrincipalRayModel_CFoV.m
    * ResDistSPC_CFoV.m
    * StitchedResDistLytro.m
    * StitchedResDistR29.m
    * ERR.m
  
`Base_m.m` is a function for projecting a set of light cones on 
a given base plane. This functionality is used e.g. inside the aperture to find 
the part of the light cone that can pass through the aperture geometry.

`ResDistPrincipalRayModel_CFoV.m` calculates the minimum resolvability distance 
in the common field of view of all the microlenses using the principal-ray-model. 
The resolvability distance is calculated for several depth planes.
The principal ray model is implemented using the SPC model, when the aperture 
is at the center of the microlenses and in the size of a pin-hole. 
The results have been used in papers [SPIE13] and [Raytrix].
More info about the parameters can be found in the header of the file.

`ResDistSPC_CFoV.m` calculates the minimum resolvability distance in the common 
field of view (CFoV) of all the microlenses using the SPC model. 
The resolvability distance is calculated for several depth planes.
More info about the parameters can be found in the header of the file.

`StitchedResDistLytro.m` generates the SPC model of a Lytro plenoptic camera and 
then calculates the minimum-resolvability-distance for several depth planes. 
The surface area to look for the  minimum-resolvability-distance goes outside 
the common field of view (CFoV), using an incremental approach. 
This is described in detail in paper [Lytro].

`StitchedResDistR29.m` generates the SPC model of an R29 plenoptic camera and 
then calculates the minimum-resolvability-distance for several depth planes. 
The surface area to look for the  minimum-resolvability-distance goes outside 
the common field of view (CFoV), using an incremental approach. 
This is described in detail in papers [Lytro] and [Raytrix].

`ERR.m`, the Effective Resolution Ratio calculates and plots the 
normalized-spatial-resolution for R29 plenoptic camera using the 
minimum-resolvability-distance values previously calculated in 
`StitchedResDistR29.m`.
It also calculates and plots the normalized-spatial-resolution for a plenoptic 
camera with three-focal length microlens array structure, using the analytical 
approach.

More info about the SPC model can be found in [a], [b] and [c].

References
-------------------------------------------------------------------------------

[SPIE13] Damghanian, Mitra, Roger Olsson, Mårten Sjöström, Hector Navarro Fructuoso, 
and Manuel Martinez-Corral. "Investigating the lateral resolution in a plenoptic 
capturing system using the SPC model." In IS&T/SPIE Electronic Imaging, pp. 
86600T-86600T. International Society for Optics and Photonics, 2013.
URL: http://miun.diva-portal.org/smash/record.jsf?pid=diva2:607155

[Raytrix] Damghanian, Mitra, Roger Olsson, Mårten Sjöström, Arne Erdmann, and 
Christian Perwass. "SPATIAL RESOLUTION IN A MULTI-FOCUS PLENOPTIC CAMERA." (2014): 
1932-1936. URL: http://miun.diva-portal.org/smash/record.jsf?pid=diva2:762088

[Lytro] Damghanian, Mitra, Roger Olsson, and Mårten Sjöström. "Performance analysis 
in Lytro camera: Empirical and model based approaches to assess refocusing quality."
ICASSP, IEEE International Conference on Acoustics, Speech and Signal Processing - 
Proceedings, IEEE conference proceedings, 2014, 559-563 (2014).
URL: http://miun.diva-portal.org/smash/record.jsf?pid=diva2:706716

[a] Damghanian, Mitra, Roger Olsson, and Mårten Sjöström. "The Sampling Pattern 
Cube–A Representation and Evaluation Tool for Optical Capturing Systems." In 
Advanced Concepts for Intelligent Vision Systems, pp. 120-131. Springer Berlin 
Heidelberg, 2012. URL: http://miun.diva-portal.org/smash/record.jsf?pid=diva2:555237

[b] Damghanian, Mitra, Roger Olsson, and Marten Sjostrom. "Extraction of the lateral 
resolution in a plenoptic camera using the SPC model." In 3D Imaging (IC3D), 2012 
International Conference on, pp. 1-5. IEEE, 2012.
URL: http://miun.diva-portal.org/smash/record.jsf?pid=diva2:582321

[c] Damghanian, Mitra. “The sampling pattern cube: A framework for representation 
and evaluation of plenoptic capturing systems,” Licentiate thesis No 99, Mid
Sweden University, Department of Information and Communication Systems, 2013.
URL: http://miun.diva-portal.org/smash/record.jsf?pid=diva2:626787
