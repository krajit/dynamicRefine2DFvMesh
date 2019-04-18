## Adaptive mesh refinement in OpenFoam-v1812 for 2-dimensional problems

Standard OpenFOAM-v1812 has a module called dynamicRefineFvMesh. This adaptive mesh refinement tool does not seem to be 
applicable for two dimensional problems. [Shonibare](https://www.academia.edu/16217705/Two-dimensional_adaptive_meshing_in_OpenFOAM) and group 
recently described how to adapt this adaptive mesh tools for two dimensional problems. However, 
CFD community struggled to implement their ideas. See this [cfd-online forum](https://www.cfd-online.com/Forums/openfoam-community-contributions/118870-2d-adaptive-mesh-refinement.html). 

Eventually, [Luca Cornolti)(https://www.cfd-online.com/Forums/openfoam-community-contributions/118870-2d-adaptive-mesh-refinement.html) shared a working for code
for adaptive mesh refinement for 2d problems. His initial code were for foam-extend. This repository is simple adaptation Luca's work
for OpenFOAM-v1812. The animated screenshot clearly shows that it works. 

[Disclaimer: It is to be noted that in the same forum discussion, Luca Cornolti also shared his adaptive mesh refinement 
for OpenFOAM-v1816. This repository is not a copy of that. Instead, this codes is built on top of Luca's foam-extend codes]


<table>
    <tr>
        <td>
          <img src="https://github.com/krajit/dynamicRefine2DFvMesh/blob/master/tutorial/interFoam/damBreak/meshRefinement-2d.gif?raw=true?raw=true" alt="Smiley face" height="200px" width="550px">
        </td>        
        <td>
          <img src="https://github.com/krajit/dynamicRefine2DFvMesh/blob/master/tutorial/interFoam/damBreak/alpha-with-refinement.gif?raw=true" alt="Smiley face" height="200px" width="550px">
        </td>        
        <td>
          <img src="https://github.com/krajit/dynamicRefine2DFvMesh/blob/master/tutorial/interFoam/damBreak/meshRefinement-3d.gif?raw=true?raw=true" alt="Smiley face" height="200px" width="550px">
        </td>        
</tr>
</table>

### Compilation

Switch to src directory and run
<pre>
wmake libso
</pre>

### Example
To run the a dam-break example, switch to tutorial/interFoam/dambreak and then run
<pre>
./Allrun
</pre>
