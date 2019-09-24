# A Tool for Kinematic Analysis on MaNGA IFS Data
Fitting kinemetry to a galaxy with a strong bar (9886-12705):
<img src="https://github.com/beckynevin/MaNGA_Kinemetry_Tool/blob/master/kinemetry_figures/kinemetry_result_9886-12705.png" width=600>


# What it does

This tool utilizes Marvin and python command-line interfaces to:
1) Download a map from a galaxy of your choice,
2) Prepare a .txt file for input into kinemetry,
3) Write idl code to run kinemetry, and
4) Uses the output of kinemetry to extract the higher order moments and best fit disk for the galaxy.

# Who it is for
The purpose of this tool is to provide a more accessible way for users of IFS data to carry out a quick kinematic analysis.

# The specifics of this repository
The main aspect of this package is 'Extract_Kinemetry_txt.ipynb' which is a python code that downloads MaNGA maps using Marvin. It allows the user to prepare a .txt file called 'kinemetry_input_plateid.txt' that has the x and y coordinates of each spaxel and the velocity, velocity error, velocity dispersion, and velocity dispersion error at that position. More specific details are within the readme file in the kinemetry_input_txt directory.

I then include a bunch of useful idl code, some of it enables the user to run kinemetry in series on multiple galaxies. 'run_kinemetry_on_multiple_galaxies.pro' is the code that goes into the specifics of kinemetry, while 'idl_wrapper_kinemetry.pro' is the wrapper that can be called to idl along with a list of plateifus and it will initialize kinemetry to run on all of those galaxies.

For now, I'm also including write_idl_input_code.py, which writes the 'run_kinemetry_on_multiple_galaxies.pro' file, but it may be easier to simply copy the .pro file.
