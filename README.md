Misorientation
==============

Misorientation is a script that allows you to determine misorientation axis and angle across a grain boundary between two crystals A and B of known orientation (2 sets of Euler angles) and crystal parameters (a,b,c, alpha, beta, gamma). 

## Requirements
Python 2.7 with Numpy, Matplotlib, Tkinter, PIL available through Enthought Canopy or Anaconda for instance. Mac Users will need Active Tcl installed. 

# User Guide


The two crystals are represented on a stereographic projection, as shown below. 

![img1](/img1.png?raw=true)

## Plot procedure
* Enter the crystal structure or import it from the menu bar. The structure can be modified/added by modifying the structure.txt file. The format is: name a b c alpha beta gamma space group. 
* Enter the maximum indice to define the number of poles/directions.
* Enter the Euler angles of crystal A and B
* Click the Plot button
 
## Misorientation procedure

* Click on the misorientation button with "Show indices" button ticked. The misorientation axes will be shown

![img2](/img2.png?raw=true)

* By ticking "Show angle", "Show axes" and "Show numbers" followed by "Misorientation" will show respectively the misorientation angle, axes indices and a number. All these data can be exported in a txt file, located in the parent directory, using the menu "Save" and "Save data".

## Additional plotting features
* Modify the number of poles/directions shown by increasing the d value (interplanar distance) by a given increment.
* Draw directions instead of poles by ticking the uvw button before plotting.
*  Add a pole/direction by clicking the "Add" button, or a family of equivalent poles/directions (circular permutations of the indices) by clicking "Symmetry"
*  Draw a plane by clicking the "Plane" button.

(For hexagonal structure, the third indice is omitted for poles).
