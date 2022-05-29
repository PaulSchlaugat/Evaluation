import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
from matplotlib.cbook import get_sample_data
from lmfit import Model
from lmfit.models import PowerLawModel
from scipy import signal
import scipy.constants as sc

plt.rcParams.update({
     'text.latex.preamble' : [r"\usepackage{siunitx}", r'\usepackage{amsmath}'],
     "text.usetex": True,
     "font.family": "sans-serif",
     "font.sans-serif": ["Computer Modern Roman"]})

def Read(data):
     '''This function reads the data from the .txt file and returns the
     wavelength, energy and intensity values in seperate arrays called after the values.'''
     
     wavelength, energy, intensity = [], [], []
     file = open(data, 'r')
     for line in file:
          this_line = line.split('\t')
          wavelength.append(float(this_line[0].replace(',', '.')))
          energy.append(float(this_line[1].replace(',', '.')))
          intensity.append(float(this_line[2].replace(',', '.').replace('\n', '')))
     return wavelength, energy, intensity

def FindClosestDatapoint(array, value):
     '''A function to find the closest datapoint (value) in a given array. The first input is the array, the second the value.'''
     
     minDist, pos = abs(array[0]-value), 0
     for i in range(len(array)):
          if abs(array[i]-value) < minDist: minDist, pos = abs(array[i]-value), i
     return pos

def Write_txt(save_place, data):
     '''A function to write the given data in a .txt file at the given folder (save_place).'''
     
     textdat = open(save_place, 'w+')
     for line in data: textdat.write(line)
     textdat.close()
     return

def PowerSeriesAllSpectra(folder_path, file_name, start=-1, end=-1, max_type='amplitude', fitting=False, cutoff=1, log=False, save=False, save_txt=False, Parameters=['', '']):
     '''This function takes power series data that is all saved in one .dat file and plots and fits it.'''
     
     energy, wavelength, allcounts, power = [], [], [], []
     file = open(folder_path + file_name, 'r')
     c = 0
     with open(folder_path + file_name) as f:
          firstline = f.readline(-1)
     for i in range(int(firstline[-3:])):
          allcounts.append([])
     for line in file:
          if c == 2:
               oneline = line[2:].split('\t')
               for i in range(len(oneline)):
                    oneline[i] = oneline[i].replace('\n', '')
                    power.append(1e6*0.06*float(oneline[i].replace(',', '.'))/7.8) #HERE FOR n AND u
          elif c > 3:
               energy.append(float(line[:11].replace(',', '.')))
               wavelength.append(float(line[12:23].replace(',', '.')))
               oneline = line[24:].split('\t')
               for i in range(len(oneline)):
                    oneline[i] = oneline[i].replace('\n', '')
                    allcounts[i].append(2*float(oneline[i].replace(',', '.'))) #times two, because of 0.5s acquisition time
          c+=1
     if start != -1: start = FindClosestDatapoint(wavelength, start)
     else: start = 0
     if end != -1: end = FindClosestDatapoint(wavelength, end)
     else: end = wavelength.index(wavelength[-1])
     energy, wavelength = energy[start:end], wavelength[start:end]
     for i in range(len(allcounts)): allcounts[i] = allcounts[i][start:end]
     
     x = wavelength
     def Func(x, amplitude, center, sigma, yvar):
               return amplitude/np.pi * sigma/((x - center)**2 + sigma**2) + yvar
     FitFunc = Model(Func, nan_policy = 'propagate', independent_vars=['x'])

     pars = FitFunc.make_params()
     pars['amplitude'].set(value=5000, min=0)
     pars['center'].set(value=772.2)
     pars['sigma'].set(value=1, min=0)
     pars['yvar'].set(value=0)
     
     amplitudes = []
     for counts in allcounts:
          result = FitFunc.fit(counts, pars, x=x)
          a, b, c, d = result.values["amplitude"], result.values["center"], result.values["sigma"], result.values["yvar"]
          yspäter=[]
          for i in x:
               yspäter.append(a/np.pi * c/((i - b)**2 + c**2) + d)
          if max_type == 'amplitude': amplitudeval = result.values['amplitude']
          elif max_type == 'counts': amplitudeval = max(yspäter)
          elif max_type == 'center': amplitudeval = result.values['center']
          amplitudes.append(amplitudeval)
     
     fig = plt.figure()
     fig.patch.set_facecolor('white')
     fig.patch.set_alpha(1)
     ax1 = fig.add_subplot(111)
     
     if log:
          ax1.set_xscale('log')
          ax1.set_yscale('log')
     
     ax1.plot(power, amplitudes, 'o')
     
     #This was specifically for our data:
     power, amplitudes = power[:-6], amplitudes[:-6]
     
     if fitting:
          x = power
          def PowerFit(Isat, x, PN): return Isat/(1 + PN/(x))
          PS = Model(PowerFit, nan_policy = 'propagate', independent_vars=['x'])
          
          pars = PS.make_params()
          pars['Isat'].set(value=20000, min=0)
          pars['PN'].set(value=10)
          result = PS.fit(amplitudes, pars, x=x)
          print(result.fit_report(min_correl=0.5))
          ax1.plot(x, result.best_fit, '--')
          
          cutoff = FindClosestDatapoint(x, cutoff)
          x = x[:cutoff]
          amplitudes = amplitudes[:cutoff]
          mod = PowerLawModel(nan_policy = 'propagate', independent_vars=['x'])
          pars = mod.guess(amplitudes, x=x)
          result2 = mod.fit(amplitudes, pars, x=x)
          print(result2.fit_report(min_correl=0.25))
          
          ax1.plot(x, result2.best_fit)
          
          ax1.set_ylabel(r'Counts (1/s)', size=14)
          ax1.set_xlabel(r'Excitation Power Density (W/cm$^2$)', size=14)
          if Parameters[0] == '': Parameters[0] = 'Saturation Curve'
          else: Parameters[0] = 'Saturation Curve\n' + Parameters[0]
          if Parameters[1] == '': Parameters[1] = 'Power Law'
          else: Parameters[1] = 'Power Law\n' + Parameters[1]
          ax1.legend(['Data', Parameters[0], Parameters[1]], fontsize=12)
     
     savefile = str(folder_path) + "imgs\\" + file_name + 'ENTIRESERIES'
     if log: savefile = savefile + 'log'
     if save: plt.savefig(savefile + '.pdf', dpi=300, bbox_inches='tight')
     if save_txt: Write_txt(savefile + '.txt', result.fit_report(min_correl=0.5) + result2.fit_report(min_correl=0.5))
     plt.show()
     return

def PowerSeries(folder_path, data, start=-1, end=-1, c1=770, c2=-1, a1=1000, yvars=600, Mod='Lorentz', save=False, save_txt=False, display=False, display_num=0, xtype='wavelength', max_type='amplitude', show=False, log=False, fitting=True, isat=2000, pn=1, PowerlawToken=False, timescale=1, entire=False, cutoff=-1, Parameters=['', '']):
     '''This function takes multiple files that make up a power series and plots and fits them.'''
     
     if display:
          PlotOne(folder_path=folder_path, file_name=data[display_num], start=start, end=end, xtype=xtype)
          return
     amplitudes, allxval, merken = [], [], []
     for dats in data:
          thisval = ''
          c = -7
          result = Fit(folder_path=folder_path, file_name=dats, PowerSeriesToken=True, Mod=Mod, save=False, save_txt=False, start=start, end=end, c1=c1, c2=c2, a1=a1, yvars=yvars, xtype=xtype, max_type=max_type, show=show, timescale=timescale, show_linewidth=False)
          
          while dats[c] != '_':
               thisval = dats[c] + thisval
               c-=1
          thisval = thisval.replace('p', '.')
          thisval = 1e3*0.06*float(thisval)/7.8
          merken.append(thisval)
          if allxval != []:
               Happened = False
               for i in range(0, len(allxval)-1, 1):
                    if thisval < allxval[i]:
                         allxval = allxval[:i] + [thisval] + allxval[i:]
                         amplitudes = amplitudes[:i] + [result] + amplitudes[i:]
                         Happened = True
                         break
               if Happened == False:
                    allxval.append(thisval)
                    amplitudes.append(result)
          else:
               allxval.append(thisval)
               amplitudes.append(result)
     if PowerlawToken: return [allxval, amplitudes]
     
     fig = plt.figure()
     fig.patch.set_facecolor('white')
     fig.patch.set_alpha(1)
     ax1 = fig.add_subplot(111)
     if log:
          ax1.set_xscale('log')
          ax1.set_yscale('log')
     
     ax1.plot(allxval, amplitudes, 'o')
     
     if fitting:
          x = allxval
          def PowerFit(I_sat, x, P_N):
               return I_sat * (1/(1 + P_N/x))
          PS = Model(PowerFit, nan_policy = 'propagate', independent_vars=['x'])
          
          pars = PS.make_params()
          pars['I_sat'].set(value=isat, min=0)
          pars['P_N'].set(value=pn)
          result = PS.fit(amplitudes, pars, x=x)
          print(result.fit_report(min_correl=0.5))
          ax1.plot(allxval, result.best_fit, '--')

          if Parameters[0] == '': Parameters[0] = 'Saturation Curve'
          else: Parameters[0] = 'Saturation Curve\n' + Parameters[0]
          if Parameters[1] == '': Parameters[1] = 'Power Law'
          else: Parameters[1] = 'Power Law\n' + Parameters[1]
          
          if entire:
               xdata = allxval[:FindClosestDatapoint(allxval, cutoff)]
               ydata = amplitudes[:FindClosestDatapoint(allxval, cutoff)]
               mod = PowerLawModel(nan_policy = 'propagate', independent_vars=['x'])
               pars = mod.guess(ydata, x=xdata)
               result2 = mod.fit(ydata, pars, x=xdata)
               print(result2.fit_report(min_correl=0.5))
               ax1.plot(xdata, result2.best_fit)
               ax1.legend(['Data', Parameters[0], Parameters[1]], fontsize=12)
          else:
               ax1.legend(['Data', Parameters[0]], fontsize=12)
          
          ax1.set_ylabel(r'Counts (1/s)', size=14)
          ax1.set_xlabel(r'Excitation Power Density (W/cm$^2$)', size=14)
     
     savefile = str(folder_path) + 'imgs\\ENTIRESERIES'
     if log: savefile = savefile + 'log'
     if save: plt.savefig(savefile + '.pdf', dpi=300, bbox_inches='tight')
     if save_txt: Write_txt(savefile + '.txt', result.fit_report(min_correl=0.5) + result2.fit_report(min_correl=0.5))
     plt.show()
     return allxval, amplitudes

#Not needed right now
def PowerLaw(folder_path, data, cutoff=-1, start=-1, end=-1, c1=770, a1=1000, yvars=600, Mod='Lorentz', save=False, display=False, display_num=0, xtype='wavelength', max_type='amplitude', show=False, log=False, xticks=[], yticks=[], timescale=1):
     values = PowerSeries(folder_path, data, start=start, end=end, c1=c1, a1=a1, yvars=yvars, Mod=Mod, display=display, display_num=display_num, xtype=xtype, max_type=max_type, show=show, PowerlawToken=True, timescale=timescale)
     x, y = values[0], values[1]
     
     fig = plt.figure()
     fig.patch.set_facecolor('white')
     fig.patch.set_alpha(1)
     ax1 = fig.add_subplot(111)
     ax1.set_ylabel(r'Counts (1/s)', size=14)
     ax1.set_xlabel(r'Power (nW)', size=14)
     #ax1.set_xlabel(r'Power (' + chr(181) + r'W)', size=12)
     
     if log:
          ax1.set_xscale('log')
          ax1.set_yscale('log')
     
     if cutoff != -1:
          cutoff = FindClosestDatapoint(x, cutoff)
          x = x[:cutoff]
          y = y[:cutoff]
     ax1.plot(x, y, 'o')
     
     mod = PowerLawModel(nan_policy = 'propagate', independent_vars=['x'])
     pars = mod.guess(y, x=x)
     result = mod.fit(y, pars, x=x)
     print(result.fit_report(min_correl=0.25))
     ax1.plot(x, result.best_fit)
     if xticks != []: ax1.set_xticks(xticks)
     if yticks != []: ax1.set_yticks(yticks)
     ax1.legend(['Data', 'Power Law'], fontsize=12)
     
     savefile = str(folder_path) + "imgs\\" + str(data[0][66:-4]) + 'POWERLAW'
     if log: savefile = savefile + 'log'
     if save: plt.savefig(savefile + '.pdf', dpi=300, bbox_inches='tight')
     plt.show()
     return

#See if this can be integrated were needed
def PlotFunc(x, y, legend, savefile, save, xdata=[], ydata=[], xtype="wavelength", vlines=[], inset=None):
     fig, ax1 = plt.subplots() #fig.add_subplot(111)
     fig.patch.set_facecolor('white')
     fig.patch.set_alpha(1)
     ax2 = ax1.twiny()
     
     if ydata != []: ax1.plot(x, y, ".", markersize=6) #If we have a fit given -> Data as points
     else: ax1.plot(x, y) #If no fit is given -> Data as line
     c = 0
     for dat in ydata:
          if c==0:
               color = 'red'
               style = 'solid'
          elif c==1:
               color = 'blue'
               style = 'dashed'
          elif c==2:
               color = 'green'
               style = 'dashdot'
          elif c==3:
               color = 'black'
               style = 'dotted'
          ax1.plot(xdata, dat, color, linestyle=style)
          c+=1
     ax1.set_ylabel(r'Counts (1/s)', size=14)
     if xtype=="wavelength": ax1.set_xlabel(r"Wavelength (nm)", size=14)
     elif xtype=="energy": ax1.set_xlabel(r"Energy (meV)", size=14)
     new_tick_locations = np.array([min(x) + 1/12 * (max(x)-min(x)), min(x) + 3/12 * (max(x)-min(x)), min(x) + 5/12 * (max(x)-min(x)), min(x) + 7/12 * (max(x)-min(x)), min(x) + 9/12 * (max(x)-min(x)), min(x) + 11/12 * (max(x)-min(x))])
     
     tick_function = lambda t: (sc.h*(sc.c)/(t*1e-9))/(sc.e)
     ax2.set_xlim(ax1.get_xlim())
     ax2.set_xticks(new_tick_locations)
     axis2 = tick_function(new_tick_locations)
     
     for i in vlines:
          ax1.vlines(i, min(y), max(y), 'red', '--', linewidth=1.3)
     
     if xtype=="wavelength":
          for i in range(len(axis2)): axis2[i] = (int(axis2[i]*10000))/10
     elif xtype=="energy":
          for i in range(len(axis2)): axis2[i] = int(axis2[i] *10000)/10
     ax2.set_xticklabels(axis2)
     if xtype=="wavelength": ax2.set_xlabel(r"Energy (meV)", size=14)
     elif xtype=="energy": ax2.set_xlabel(r"Wavelength (nm)", size=14)
     ax1.legend(legend, fontsize=12, loc='upper right')
     
     if inset != None:
          fn = get_sample_data(inset[0], asfileobj=False)
          arr_img = plt.imread(fn, format='png')

          imagebox = OffsetImage(arr_img, zoom=inset[2])
          imagebox.image.axes = ax1

          ab = AnnotationBbox(imagebox, inset[1])
          ax1.add_artist(ab)
     
     if save: plt.savefig(savefile + ".pdf", dpi=300, bbox_inches='tight')
     plt.show()
     return

def PlotOne(folder_path, file_name, save=False, start=0, end=10000, xtype="wavelength", vlines=[], timescale=1, legend=[], inset=None):
     '''This function plots the data given in one file.'''
     #Load the data:
     Data=Read(data=file_name)
     wavelength, energy, intensity = Data[1], Data[0], Data[2]
     for i in range(len(energy)):
          energy[i] = 1000*energy[i]
          intensity[i] = timescale*intensity[i]
     
     #Define start and end points:
     if xtype == 'wavelength':
          start = FindClosestDatapoint(wavelength, start)
          end = FindClosestDatapoint(wavelength, end)
          wavelength, energy, intensity = wavelength[start:end], energy[start:end], intensity[start:end]
     elif xtype == 'energy':
          start = FindClosestDatapoint(energy, start)
          end = FindClosestDatapoint(energy, end)
          wavelength, energy, intensity = wavelength[end:start], energy[end:start], intensity[end:start]
     
     #Create saving-path:
     savefile = str(folder_path) + "imgs\\" + "Plot" + str(file_name[0:len(file_name)-4])
     savefile = savefile + "s" + str(start) + "e" + str(end)
          
     #Plot the data (with intensity over wavelength aswell as energy):
     if legend == []: legend = ["Data"]
     if xtype=="wavelength": x, x2 = wavelength, energy
     else: x, x2 = energy, wavelength
     PlotFunc(x, intensity, legend, savefile, save, xtype=xtype, vlines=vlines, inset=inset)
     
     #Printing the Name:
     Name, counter = '', 0
     for letter in file_name:
          if counter > 6: Name = Name + letter
          if letter == "\\": counter+=1
     if counter == 0: print(file_name)
     else: print(Name)
     #print("Peak position: " + str(xval) + Einheit + str(ymax) + "Int.")
     #if Int: return ymax-min(intensity)
     return

def Fit(file_name, folder_path, Mod, save=False, save_txt=False, start=-1, end=-1, c1=770, c2=-1, c3=770, a1=1000, a2=500, mini=0, yvars=600, xtype="wavelength", PowerSeriesToken=False, max_type='amplitude', double=False, show=True, show_linewidth=True, timescale=1, Parameters='', Temp_Series=False):
     '''With this function peaks in a spectrum can be fitted with different peak shapes.'''
     
     Data=Read(file_name)
     wavelength, energy, intensity = Data[1], Data[0], Data[2]
     for i in range(len(intensity)):
          intensity[i] = timescale*intensity[i]
     
     for i in range(len(energy)):
          energy[i] = 1000*energy[i]
     
     if xtype == 'wavelength':
          start = FindClosestDatapoint(wavelength, start)
          end = FindClosestDatapoint(wavelength, end)
          wavelength, energy, intensity = wavelength[start:end], energy[start:end], intensity[start:end]
     elif xtype == 'energy':
          start = FindClosestDatapoint(energy, start)
          end = FindClosestDatapoint(energy, end)
          wavelength, energy, intensity = wavelength[end:start], energy[end:start], intensity[end:start]
     
          
     savefile = str(folder_path) + "imgs\\" + "Fit" + str(file_name[0:len(file_name)-4]) + Mod
     savefile = savefile + "s" + str(start) + "e" + str(end)
     
     x = wavelength
     xenergy = energy
     linewidth1 = 0
     
     if c1==-1: c1=x[intensity.index(max(intensity))]
     if c2==-1: c2=c1+1
     wavelength, energy, intensity = np.asarray(wavelength), np.asarray(energy), np.asarray(intensity)
     xspäter=np.linspace(min(x), max(x), 1000)
     xenergyspäter=np.linspace(max(xenergy), min(xenergy), 1000)
     
     if Mod=="Lorentz with Gauss":
          def Func(x, amplitude1, center1, sigma1, amplitude2, center2, sigma2, yvar):
               return amplitude1/np.pi * sigma1/((x - center1)**2 + sigma1**2) + amplitude2/(sigma2*np.sqrt(2*np.pi)) * np.e**((-(x-center2)**2)/(2*sigma2**2)) + yvar

          FitFunc = Model(Func, nan_policy = 'propagate', independent_vars=['x'])

          pars = FitFunc.make_params()
          pars['amplitude1'].set(value=a1, min=1)
          pars['center1'].set(value=c1, min=c1-1, max=c1+1)
          pars['sigma1'].set(value=1, min=0)
          pars['amplitude2'].set(value=a2, min=0)
          pars['center2'].set(value=c2, min=c2-1, max=c2+1)
          pars['sigma2'].set(value=1, min=0)
          pars['yvar'].set(value=yvars, min=mini)
          result = FitFunc.fit(intensity, pars, x=x)
          if show and not Temp_Series: print(result.fit_report(min_correl=1)) #0.25
          
          a, b, c, d, e, f, g = result.values["amplitude1"], result.values["center1"], result.values["sigma1"], result.values["amplitude2"], result.values["center2"], result.values["sigma2"], result.values["yvar"]
          
          l1, l2 = result.values["center1"] - result.values["sigma1"], result.values["center1"] + result.values["sigma1"]
          linewidth1 = ((sc.h * sc.c * (1/(l1*10**(-9)) - 1/(l2*10**(-9))))*6.242*10**(18))*1e6
          # l1, l2 = result.values["center2"] - result.values["sigma2"], result.values["center2"] + result.values["sigma2"]
          # linewidth2 = ((sc.h * sc.c * (1/(l1*10**(-9)) - 1/(l2*10**(-9))))*6.242*10**(18))*1e6
          
          yspäter=[]
          yspäterl=[]
          yspäterg=[]
          for i in xspäter:
               yspäter.append(a/np.pi * c/((i - b)**2 + c**2) + d/(f*np.sqrt(2*np.pi))*np.e**((-(i-e)**2)/(2*f**2)) + g)
               yspäterl.append(a/np.pi * c/((i - b)**2 + c**2) + g)
               yspäterg.append(d/(f*np.sqrt(2*np.pi)) * np.e**((-(i-e)**2)/(2*f**2)) + g)
               
          ydata=[yspäter, yspäterl, yspäterg]
          if show_linewidth and not Temp_Series: print('Linewidth: ' + str(int(linewidth1))) #legend=["Data", "Fit", "Lorentzian " + str(int(linewidth1)) + "ueV", "Gaussian"]
          
          if Parameters == '': Parameters = 'Lorentzian'
          else: Parameters = 'Lorentzian\n' + Parameters
          
          legend=["Data", "Fit", Parameters, "Gaussian"]
          
          if xtype=="wavelength": x = wavelength #THIIIIISSSS
          elif xtype=="energy":
               x = energy
               xspäter = xenergyspäter
          if show: PlotFunc(x, intensity, legend, savefile, save, xdata=xspäter, ydata=ydata, xtype=xtype)
          if Temp_Series: return [x, intensity, xspäter, yspäter, linewidth1, file_name, result.values['center1'], result.values['sigma1']]
          
          if double:
               if max_type == 'amplitude': amplitudeval = [result.values['amplitude1'], result.values['amplitude2']]
               elif max_type == 'counts': amplitudeval = [max(yspäterl), max(yspäterg)]
               elif max_type == 'center': amplitudeval = [result.values['center1'], result.values['center2']]
          else:
               if max_type == 'amplitude': amplitudeval = result.values['amplitude1']
               elif max_type == 'counts': amplitudeval = max(yspäter)
               elif max_type == 'center': amplitudeval = result.values['center1']
               
          Name, counter = '', 0
          for letter in file_name:
               if counter > 6: Name = Name + letter
               if letter == "\\": counter+=1
          if counter == 0: print(file_name)
          else: print(Name)
          
          if PowerSeriesToken: return amplitudeval
     elif Mod=="Lorentz":
          def Func(x, amplitude, center, sigma, yvar):
               return amplitude/np.pi * sigma/((x - center)**2 + sigma**2) + yvar

          FitFunc = Model(Func, nan_policy = 'propagate', independent_vars=['x'])

          pars = FitFunc.make_params()
          pars['amplitude'].set(value=a1, min=1)
          pars['center'].set(value=c1)
          pars['sigma'].set(value=1, min=0)
          pars['yvar'].set(value=yvars, min=mini)
          result = FitFunc.fit(intensity, pars, x=x)
          if show and not Temp_Series: print(result.fit_report(min_correl=0.5))
          
          a, b, c, d = result.values["amplitude"], result.values["center"], result.values["sigma"], result.values["yvar"]
          l1, l2 = result.values["center"] - result.values["sigma"], result.values["center"] + result.values["sigma"]
          linewidth1 = ((sc.h * sc.c * (1/(l1*10**(-9)) - 1/(l2*10**(-9))))*6.242*10**(18))*1e6
          
          yspäter=[]
          for i in xspäter:
               yspäter.append(a/np.pi * c/((i - b)**2 + c**2) + d)
          
          ydata=[yspäter]
          
          if Parameters == '': Parameters = 'Lorentzian'
          else: Parameters = 'Lorentzian\n' + Parameters
          
          if show_linewidth and not Temp_Series: print(int(linewidth1))
          
          if Parameters == '': Parameters = 'Lorentzian'
          else: Parameters = 'Lorentzian\n' + Parameters
          
          legend = ["Data", Parameters]
          
          if xtype=="wavelength": x = wavelength #THIIIIISSSS
          elif xtype=="energy":
               x = energy
               xspäter = xenergyspäter
          
          if show: PlotFunc(x, intensity, legend, savefile, save, xspäter, ydata)
          if Temp_Series: return [x, intensity, xspäter, yspäter, linewidth1, file_name, result.values['center'], result.values['sigma']]
          
          if max_type == 'amplitude': amplitudeval = result.values['amplitude']
          elif max_type == 'counts': amplitudeval = max(yspäter)
          elif max_type == 'center': amplitudeval = result.values['center']
          
          Name, counter = '', 0
          for letter in file_name:
               if counter > 6: Name = Name + letter
               if letter == "\\": counter+=1
          #if counter == 0: print(file_name)
          #else: print(Name)
          if show:
               if counter == 0: print(file_name)
               else: print(Name)
          if show: print("Linewidth: " + str(linewidth1) + "ueV or " + str(2*float(result.values["sigma"])) + "nm at " + str(result.values["center"]) + "nm")
     elif Mod=="Gauss":
          def Func(x, amplitude, center, sigma, yvar):
               return amplitude/(sigma * np.sqrt(2 * np.pi)) * np.e**((-(x-center)**2)/(2*sigma**2)) + yvar

          FitFunc = Model(Func, nan_policy = 'propagate', independent_vars=['x'])

          pars = FitFunc.make_params()
          pars['amplitude'].set(value=a1, min=1)
          pars['center'].set(value=c1)
          pars['sigma'].set(value=1, min=0)
          pars['yvar'].set(value=yvars, min=mini)
          result = FitFunc.fit(intensity, pars, x=x)
          if show: print(result.fit_report(min_correl=0.5))
          
          a, b, c, d = result.values["amplitude"], result.values["center"], result.values["sigma"], result.values["yvar"]
          l1, l2 = result.values["center"] - result.values["sigma"], result.values["center"] + result.values["sigma"]
          linewidth1 = ((sc.h * sc.c * (1/(l1*10**(-9)) - 1/(l2*10**(-9))))*6.242*10**(18))*1e6
          
          yspäter=[]
          for i in xspäter:
               yspäter.append(a/np.pi * c/((i - b)**2 + c**2) + d)
          
          ydata=[yspäter, yspäter]
          
          if Parameters == '': Parameters = 'Gaussian'
          else: Parameters = 'Gaussian\n' + Parameters
          
          legend = ["Data", Parameters]
          
          if Temp_Series: return [x, intensity, xspäter, yspäter, linewidth1, file_name]
          if show: PlotFunc(x, intensity, legend, savefile, save, xspäter, ydata)
          
          Name, counter = '', 0
          for letter in file_name:
               if counter > 6: Name = Name + letter
               if letter == "\\": counter+=1
          if counter == 0: print(file_name)
          else: print(Name)
          if show: print("Linewidth: " + str(linewidth1) + "ueV or " + str(2*float(result.values["sigma"])) + "nm at " + str(result.values["center"]) + "nm")
     elif Mod=="Double Lorentz with Gauss":
          def Func(x, amplitude1, center1, sigma1, amplitude2, center2, sigma2, amplitude3, center3, sigma3, yvar):
               return amplitude1/np.pi * sigma1/((x - center1)**2 + sigma1**2) + amplitude2/np.pi * sigma2/((x - center2)**2 + sigma2**2) + amplitude3/(sigma3 * np.sqrt(2 * np.pi)) * np.e**((-(x-center3)**2)/(2*sigma3**2)) + yvar

          FitFunc = Model(Func, nan_policy = 'propagate', independent_vars=['x'])

          pars = FitFunc.make_params()
          pars['amplitude1'].set(value=a1, min=1)
          pars['center1'].set(value=c1)
          pars['sigma1'].set(value=0.1, min=0)
          pars['amplitude2'].set(value=a2, min=1)
          pars['center2'].set(value=c2)
          pars['sigma2'].set(value=0.4, min=0)
          pars['amplitude3'].set(value=50, min=1, max=500)
          pars['center3'].set(value=c3)
          pars['sigma3'].set(value=0.1, min=0)
          pars['yvar'].set(value=yvars, min=mini)
          result = FitFunc.fit(intensity, pars, x=x)
          if show and not Temp_Series: print(result.fit_report(min_correl=0.5))
          
          a, b, c, d, e, f, g, h, k, j = result.values["amplitude1"], result.values["center1"], result.values["sigma1"], result.values["amplitude2"], result.values["center2"], result.values["sigma2"], result.values["amplitude3"], result.values["center3"], result.values["sigma3"], result.values["yvar"]
          l1, l2 = result.values["center1"] - result.values["sigma1"], result.values["center1"] + result.values["sigma1"]
          linewidth1 = ((sc.h * sc.c * (1/(l1*10**(-9)) - 1/(l2*10**(-9))))*6.242*10**(18))*1e6
          # l1, l2 = result.values["center2"] - result.values["sigma2"], result.values["center2"] + result.values["sigma2"]
          # linewidth2 = ((sc.h * sc.c * (1/(l1*10**(-9)) - 1/(l2*10**(-9))))*6.242*10**(18))*1e6
          
          yspäter=[]
          yspäter1=[]
          yspäter2=[]
          yspäter3=[]
          for i in xspäter:
               yspäter.append(a/np.pi * c/((i - b)**2 + c**2) + d/np.pi * f/((i - e)**2 + f**2) + g/(k*np.sqrt(2*np.pi)) * np.e**((-(i-h)**2)/(2*k**2)) + j)
               yspäter1.append(a/np.pi * c/((i - b)**2 + c**2) + j)
               yspäter2.append(d/np.pi * f/((i - e)**2 + f**2) + j)
               yspäter3.append(g/(k*np.sqrt(2*np.pi)) * np.e**((-(i-h)**2)/(2*k**2)) + j)
          
          ydata=[yspäter, yspäter1, yspäter2, yspäter3]
          if show_linewidth: print(int(linewidth1*10)/10)
          
          if Parameters == '': Parameters = 'Lorentzian 1'
          else: Parameters = 'Lorentzian 1\n' + Parameters
          
          legend = ["Data", "Fit", Parameters, "Lorentzian 2", "Gaussian"]
          
          if xtype=="wavelength": x = wavelength #THIIIIISSSS
          elif xtype=="energy":
               x = energy
               xspäter = xenergyspäter
          if show: PlotFunc(x, intensity, legend, savefile, save, xspäter, ydata, xtype=xtype)
          if Temp_Series: return [x, intensity, xspäter, yspäter, linewidth1, file_name, result.values['center1'], result.values['sigma1']]
          
          Name, counter = '', 0
          for letter in file_name:
               if counter > 6: Name = Name + letter
               if letter == "\\": counter+=1
          if counter == 0: print(file_name)
          else: print(Name)
     else: raise Exception('This Model does not exist (yet)!')
     if save_txt: Write_txt(savefile + '.txt', result.fit_report(min_correl=0.5) + '\n\nLinewidth: ' + str(linewidth1) + 'ueV')
     if PowerSeriesToken: return amplitudeval
     return savefile

def Temp(folder_path, files, start=-1, end=10000, save=False, legend=['', '', '', ''], xtype='wavelength', Mod='Lorentz with Gauss', show_linewidth=False, c1=-1, c2=776, show=False, timescale=1):
     '''With this function the temperature series collected in multiple files can be plotted with the lifetimes in the legend.'''
     
     all_fits = []
     for file in files:
          print(file)
          if '_12K' in file: stats = Fit(folder_path=folder_path, file_name=file, c1=781.6, c2=780.5, c3=782.8, start=start, end=end, xtype=xtype, Mod='Double Lorentz with Gauss', show=show, Temp_Series=True, show_linewidth=show_linewidth, timescale=timescale)
          else: stats = Fit(folder_path=folder_path, file_name=file, c1=c1, c2=c2, start=start, end=end, xtype=xtype, Mod=Mod, show=show, Temp_Series=True, show_linewidth=show_linewidth, timescale=timescale)
          all_fits.append(stats)
     
     for i in all_fits:
          print(i[4], i[5], i[6], i[7])
     
     fig = plt.figure()
     fig.patch.set_facecolor('white')
     fig.patch.set_alpha(1)
     ax1 = fig.add_subplot(111)
     ax2 = ax1.twiny()
     
     ax1.plot(all_fits[0][0], all_fits[0][1])
     ax1.plot(all_fits[1][0], all_fits[1][1])
     ax1.plot(all_fits[2][0], all_fits[2][1])
     ax1.plot(all_fits[3][0], all_fits[3][1])
     ax1.legend(['Data' + legend[0], 'Data' + legend[1], 'Data' + legend[2], 'Data' + legend[3]], fontsize=12) #[p1, p2, p3, p4], 
     ax1.set_xlabel(r'Energy (meV)', size=14)
     ax1.set_ylabel(r'Counts (1/s)', size=14)
     
     x = stats[0]
     new_tick_locations = np.array([min(x) + 1/12 * (max(x)-min(x)), min(x) + 3/12 * (max(x)-min(x)), min(x) + 5/12 * (max(x)-min(x)), min(x) + 7/12 * (max(x)-min(x)), min(x) + 9/12 * (max(x)-min(x)), min(x) + 11/12 * (max(x)-min(x))])
     tick_function = lambda t: (sc.h*(sc.c)/(t*1e-9))/(sc.e)
     ax2.set_xlim(ax1.get_xlim())
     ax2.set_xticks(new_tick_locations)
     axis2 = tick_function(new_tick_locations)
     
     if xtype=="wavelength":
          for i in range(len(axis2)): axis2[i] = (int(axis2[i]*10000))/10
     elif xtype=="energy":
          for i in range(len(axis2)): axis2[i] = int(axis2[i] *10000)/10
     ax2.set_xticklabels(axis2)
     if xtype=="wavelength": ax2.set_xlabel(r"Energy (meV)", size=14)
     elif xtype=="energy": ax2.set_xlabel(r"Wavelength (nm)", size=14)
     
     savefile = folder_path + 'imgs\\Temp_Series'
     print(savefile)
     if save: plt.savefig(savefile + ".pdf", dpi=300, bbox_inches='tight')
     plt.show()
     return

def HBT(folder_path, file_name, start=0, end=0, save=False, mode='same', norms=False, offs=False, save_txt=False, binning=1, Parameters=['', ''], Jitter=0.03):
     '''In this function data from autocorrelation measurements can be fitted. Needed are the folder path to the data 
     and the file_name in the folder. With variables start and end the data can be cut at specific times. These times
     should be given as the times in the plot and not as indices in an array. Variables save and save_txt can the plot
     and fit parameters be saved. The variable mode decides about the mode of the convolution (either full or same).
     Variables norms and offs precalculate the normalization and offset from zero of the data. In the normalization the first
     and last fourth of the data is used. The offset places the minimum of the data at the time zero. The variable binning can
     be used to bin the data.'''
     
     #Jitter of respective detector:
     if Jitter == 'APD': Jitter = np.sqrt(2)*0.35
     elif Jitter == 'APDs': Jitter = 0.4
     elif Jitter == 'NW': Jitter = np.sqrt(2)*0.03
     elif Jitter == 'NWs': Jitter = 0.04
     else: Jitter = 0.05
     
     #Read the Data from file:
     time, counts = [], []
     file = open(folder_path + file_name, 'r')
     for line in file:
          if line[0] != '#':
               time.append(1e9*float(line[:15].replace(',', '.')))
               counts.append(float(line[19:33].replace(',', '.')))
     
     #Cut the data if needed
     if start == 0: start=time[0]
     if end == 0: end=time[-1]
     x = time[FindClosestDatapoint(time, start):FindClosestDatapoint(time, end)]
     y = counts[FindClosestDatapoint(time, start):FindClosestDatapoint(time, end)]
     x = np.asarray(x)
     y = np.asarray(y)
     
     #Offset:
     if offs:
          Offs = x[y.index(min(y))]
          x -= Offs
     
     #Normalization:
     if norms:
          normalization = (np.sum(y[:int(len(y)/4)]) + np.sum(y[3*int(len(y)/4):])) / (len(y[:int(len(y)/4)]) + len(y[3*int(len(y)/4):]))
          y = y/normalization
     
     #Binning:
     if binning != 1:
          #Cut the array so that each bin is same size:
          binning_overflow = -1*(len(x)%binning)
          if binning_overflow != 0: x = x[:binning_overflow]
          
          x_binned, y_binned = [], []
          for i in range(len(x)):
               if (i+1)%binning == 0:
                    x_binned.append(np.sum(x[(i+1)-binning:(i+1)])/binning)
                    y_binned.append(np.sum(y[(i+1)-binning:(i+1)])/binning)
          x, y = x_binned, y_binned
     
     #Create the output figure:
     fig = plt.figure()
     fig.patch.set_facecolor('white')
     fig.patch.set_alpha(1)
     ax1 = fig.add_subplot(111)
     
     #Plot the raw Data:
     ax1.plot(x, y)
     
     #Define the fit-function:
     def HBT_Fit(x, Norm, g2, tau, Offset):
          s2l2 = 2.3548200450309493 # = 2*np.sqrt(2*np.log(2))
          fwhm = Jitter
          sigma = fwhm/s2l2
          
          #Calculate function values:
          Af = Norm*(1-((1-g2)*np.e**(-abs(x+Offset)/tau)))
          Ag = 1/(sigma*2.5066282746310002) * np.e**(-(x)**2/(2*sigma**2))
          Ag = Ag/sum(Ag)
          
          #Wrap Af around:
          wrap = 10
          Af = np.append(np.append(Af[-int(len(Af)/wrap):], Af), Af[:int(len(Af)/wrap)])
               
          #Convolution:
          Afg = signal.fftconvolve(Af, Ag, mode=mode)
          
          #Truncating start and end if mode 'full' has been used:
          x_len = len(x)
          front = math.floor(x_len/2) + int(x_len/wrap)
          back = x_len + math.floor(x_len/2) + int(x_len/wrap)
          if mode == 'full': Afg = Afg[front:back]
          
          return Afg
     
     #Create a Model-Object of the fit function:
     FitFunc = Model(HBT_Fit, nan_policy = 'propagate', independent_vars=['x'])
     pars = FitFunc.make_params()
     pars['Norm'].set(value=1, min=0.95, max=1.05)
     pars['g2'].set(value=0.2, min=0.004, max=0.5)
     pars['tau'].set(value=0.1, min=0.01, max=15)
     pars['Offset'].set(value=0, min=-5, max=5)
     
     #Fitting:
     result = FitFunc.fit(y, pars, x=x)
     print(result.fit_report(min_correl=0.3))
     
     #Plot the fit:
     ax1.plot(x, result.best_fit)
     
     #Set labels and legend:
     ax1.set_xlabel(r'Time (ns)', size=14)
     ax1.set_ylabel(r'Normalized Counts (arb. units)', size=14)
     if Parameters[0] == '': Parameters[0] = 'Autocorrelation Fit'
     else: Parameters[0] = 'Autocorrelation Fit\n' + Parameters[0]
     ax1.legend(['Data', Parameters[0]], fontsize=12)
     
     ax1.set_xlim(-50, 50)
     ax1.set_ylim(0, 1.5)
     
     #Saving the picture:
     savefile = str(folder_path) + "\\imgs\\" + str(file_name) + 's' + str(start) + 'e' + str(end) + 'binning is ' + str(binning)
     if save: plt.savefig(savefile + '.pdf', dpi=300, bbox_inches='tight', transparent='true')
     if save_txt: Write_txt(savefile + '.txt', result.fit_report(min_correl=0.25))
     plt.show()
     print('min of y: ', min(y), '; min of fit: ', min(result.best_fit))
     return

def Lifetime(folder_path, file_name, exp_start=-1, start=-1, end=-1, exponents=2, display=False, save=False, log=False, yticks=[], txt_save=False, Besonders=False, Parameters=''):
     '''This function takes a lifetime measurement and plots and fits the data.'''
     
     #Read Data:
     time, counts = [], []
     file = open(folder_path + file_name, 'r')
     if Besonders:
          for line in file:
               this_line = line.split('\t')
               time.append(float(this_line[0].replace(',', '.')))
               counts.append(float(this_line[1].replace(',', '.')))
     else:
          for line in file:
               if line[0]!='#':
                    time.append(float(line[1:15].replace(',', '.')))
                    counts.append(float(line[19:33].replace(',', '.')))
     #Define start and end points:
     if start != -1: start = FindClosestDatapoint(time, start)
     else: start = FindClosestDatapoint(time, 0)
     if end != -1: end = FindClosestDatapoint(time, end)
     else: end = FindClosestDatapoint(time, time[-1])
     time = time[start:end]
     counts = counts[start:end]
     offset = time[0]
     for i in range(len(time)):
          time[i] = time[i] - offset
     
     for i in range(len(time)):
         time[i] = time[i]*1e9
     fig = plt.figure()
     fig.patch.set_facecolor('white')
     fig.patch.set_alpha(1)
     ax1 = fig.add_subplot(111)
     if log:
          ax1.set_yscale('log')
     ax1.plot(time, counts, 'o', markersize=3)
     for i in range(len(time)):
         time[i] = time[i]*1e-9
     #Fitting:
     if display == False:
          #Shift the time to be at 0 when the exp should start and y down for better log
          if exp_start==-1: timeoffset = time[counts.index(max(counts))]
          else: timeoffset = time[FindClosestDatapoint(time, exp_start)]
          for i in range(len(time)):
               time[i] = time[i] - timeoffset
          cut = time.index(0)
          x = time[cut:]
          counts = counts[cut:]
          # exp_model = Model(exp_model1, nan_policy = 'propagate', independent_vars=['x']) with prefixes for possible multiexponential
          def exp_model1(x, amplitude1, decay1): #, yval): #, yval
               return amplitude1*np.e**(-x*1e9/decay1) # + yval
          def exp_model2(x, amplitude1, decay1, amplitude2, decay2, yval):
               return amplitude1*np.e**(-x*1e9/decay1) + amplitude2*np.e**(-x*1e9/decay2) + yval
          
          if exponents==1: Func = exp_model1
          elif exponents==2: Func = exp_model2
          
          FitFunc = Model(Func, nan_policy = 'propagate', independent_vars=['x'])
          #Hiermit könnte ich auch eine eigene Exp Func machen mit Prefix
          
          pars = FitFunc.make_params()
          pars['amplitude1'].set(value=1000, min=0)
          pars['decay1'].set(value=5, min=0)
          #pars['yval'].set(value=min(counts), min=min(counts)-4*min(counts), max=max(counts)+4*min(counts))
          if exponents==2: pars['amplitude2'].set(value=500, min=0)
          if exponents==2: pars['decay2'].set(value=40, min=0)
          result = FitFunc.fit(counts, pars, x=x)
          print(result.fit_report(min_correl=0.5))
          
          for i in range(len(x)):
               x[i] = x[i] + timeoffset
          for i in range(len(x)):
               x[i] = x[i]*1e9
          ax1.plot(x, result.best_fit, linewidth=2.5)
          for i in range(len(x)):
               x[i] = x[i]*1e-9
          ax1.set_xticks([0, 2, 4, 6, 8, 10, 12])
          if yticks != [] and log: ax1.set_yticks(yticks)
          ax1.set_xlabel(r"Time (ns)", size=14)
          ax1.set_ylabel(r"Counts (1/s)", size=14)
          
          if Parameters == '': Parameters = 'Exponential Decay'
          else: Parameters = r'Exponential Decay' + '\n' + Parameters
          ax1.legend([r'Data', Parameters], fontsize=12)
     
     savefile = str(folder_path) + "imgs\\" + file_name + 's' + str(start) + 'e' + str(end)
     if log: savefile = savefile + '_log'
     savefile = savefile + '_Lifetime'
     if save: plt.savefig(savefile + '.pdf', dpi=300, bbox_inches='tight')
     if txt_save: Write_txt(savefile + '.txt', result.fit_report(min_correl=0.5))
     plt.show()
     return