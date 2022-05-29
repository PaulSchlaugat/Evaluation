# How to use the functions defined in the python program
This file will go through the individual functions defined in the "Evaluation.py" file, which was created to plot and fit the data of my bachelors thesis. I will roughly explain what each function does and can be used for. Additionally, examples on how to run the functions and inputs to each will be given.
The libraries imported at the beginning will be used to plot, calculate and fit the data. The documentations of all the modules can be found here:<br/>
NumPy: https://numpy.org/doc/stable/<br/>
math: https://docs.python.org/3/library/math.html<br/>
matplotlib.pyplot: https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.html<br/>
lmfit: https://lmfit.github.io/lmfit-py/model.html<br/>
scipy: https://scipy.org/<br/>

## Basic Functions:<br/>
The following three functions will be used to read, work with and save the data. They will not directly be used to plot or fit data.<br/>
1) Read(data):<br/>
This function, as the name suggests, reads an input file containing the wavelength, energy and intensity of a photoluminescence measurement. The input should be a .txt file with 3 columns and any amount of rows of data. The first column will be interpreted as the wavelength, the second as the energy and the third as the intensity of the spectrum. Each of the three is returned in a seperate array.<br/>
2) FindClosestDatapoint(array, value):<br/>
With this function the position of the closest point in an array to a chosen value is returned. This will be used to cut the spectra or other data at desired values.<br/>
3) Write_txt(save_place, data):<br/>
This function saves texts (here used to save the fit parameters) to .txt files in a desired folder. The folder has to be specified as the save_place variable.<br/>
