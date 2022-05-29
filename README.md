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

## Photoluminescence Setup:<br/>
Text.<br/>
1) PowerSeries(folder_path, data, start=-1, end=-1, c1=770, c2=-1, a1=1000, yvars=600, Mod='Lorentz', save=False, save_txt=False, display=False, display_num=0, xtype='wavelength', max_type='amplitude', show=False, log=False, fitting=True, isat=2000, pn=1, PowerlawToken=False, timescale=1, entire=False, cutoff=-1, Parameters=['', '']): <br/>
Text.<br/>
2) PowerSeriesAllSpectra(folder_path, file_name, start=-1, end=-1, max_type='amplitude', fitting=False, cutoff=1, log=False, save=False, save_txt=False, Parameters=['', '']):<br/>
Text.<br/>
3) PlotOne(folder_path, file_name, save=False, start=0, end=10000, xtype="wavelength", vlines=[], timescale=1, legend=[], inset=None):<br/>
Text.<br/>
4) Fit(file_name, folder_path, Mod, save=False, save_txt=False, start=-1, end=-1, c1=770, c2=-1, c3=770, a1=1000, a2=500, mini=0, yvars=600, xtype="wavelength", PowerSeriesToken=False, max_type='amplitude', double=False, show=True, show_linewidth=True, timescale=1, Parameters='', Temp_Series=False):<br/>
Text.<br/>
5) Temp(folder_path, files, start=-1, end=10000, save=False, legend=['', '', '', ''], xtype='wavelength', Mod='Lorentz with Gauss', show_linewidth=False, c1=-1, c2=776, show=False, timescale=1):<br/>
Text.<br/>

## Time Resolved Setup:<br/>
Text.<br/>
1) HBT(folder_path, file_name, start=0, end=0, save=False, mode='same', norms=False, offs=False, save_txt=False, binning=1, Parameters=['', ''], Jitter=0.03):<br/>
Text.<br/>
2) def Lifetime(folder_path, file_name, exp_start=-1, start=-1, end=-1, exponents=2, display=False, save=False, log=False, yticks=[], txt_save=False, Besonders=False, Parameters=''):<br/>
Text.
