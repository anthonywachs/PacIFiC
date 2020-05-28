# Usage of the post-processing tools

## DEM

### free_surface.py

<u>Location</u>: `postProcessing/DEM/free_surface.py`

<u>Description:</u> this program plots the free surface and returns the run-out distance of a granular dam-break.

<u>Requirements</u>: Paraview5.8 for pvpython, as well as the python modules skimage, scipy, os, sys, matplotlib, numpy.

<u>Dependancies</u>: uses the python files `postProcessing/DEM/cmdLineArguments.py` and `postProcessing/DEM/generate_first_last_frames.py`

<u>Usage</u>: `python $PPHOME/DEM/free_surface.py <simulation_folder1> -<option1>=<value> <simulation_folder2> -<option1>=<value> -<option2>=<value>`

<u>Options</u>:
* `-l`: specifies the label for the previous simulation folder
* `-c`: specifies the line color for the previous simulation folder
* `-lw`: specifies the line width for the previous simulation folder
* `-ls`: specifies the line style for the previous simulation folder
* `-t`: specifies the title(s) of the plot(s)
* `-run_pv`: explicits if paraview is run. Possible values: `yes` or `no`. Default is `yes`.
* `-width`: specifies the width of the rectangular domain (in cm). Default is `20`.
* `-runout_intensity`: specifies the intensity threshold that defines the run-out distance, in percentage of the maximum achievable intensity. In a top view of the final snapshot, the intensity decreases in the direction of the avalanche. The runout distance is the first time the intensity is as low ad the value of `runout_intensity`. Default is `0.05`.

<u>Example</u>:
Typing `python3.7 $PP_HOME/DEM/free_surface.py $PP_HOME/test-cases/free_surface_testCase1 -c=blue -ls=: -l=high\ aspect\ ratio $PP_HOME/test-cases/free_surface_testCase2 -ls=black -l=low\ aspect\ ratio -width=20 -run_pv=yes -runout_intensity=0.05` leads to:

![alt text](./image/surface_plots_1.png "Plot of the free surfaces")
![alt text](./image/surface_plots_2.png "Plot of the initial and final side views for case2")
![alt text](./image/surface_plots_3.png "Plot of the final top view for case2, with intensity and run-out distances")
![alt text](./image/surface_plots_4.png "Plot of the initial and final side views for case1")
![alt text](./image/surface_plots_5.png "Plot of the final top view for case1, with intensity and run-out distances")

In the top views of the final state of the avalanche, the solid line represents the vertical-integrated intensity and the dashed vertical line reprensents the run-out distance.
