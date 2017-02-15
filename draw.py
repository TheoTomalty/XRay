import data as d
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

count_axis = "Photon Count"
angle_axis = "2$\Theta$ ($^\circ$)"
composition_axis = "Composition (%Cu)"
lattice_axis = "Lattice Parameter ($\AA$)"

def error(val):
    string = str(val)
    for letter in string:
        if letter != '0' and letter != '.':
            return letter
    
    return '0'

def dec_pos(string):
    count = 0
    for letter in string:
        if letter == ".":
            return count
        count += 1

def value(val, err):
    string1 = str(val)
    s = str(err)
    string2 = (format(err, "." + s.split("e")[-1][1:]+"f") if "e" in s else str(err))
    
    decimal1 = dec_pos(string1)
    decimal2 = dec_pos(string2)
    
    first_nonzero = 0
    for letter in string2:
        if letter != '0' and letter != '.':
            break
        first_nonzero += 1
    
    if first_nonzero < decimal2:
        index = first_nonzero - decimal2 + 1
    else:
        index = first_nonzero - decimal2 + 1
    
    return string1[:decimal1 + index] + '(' + string2[first_nonzero] + ')'

print value(13.7, 1.4)

def set_residual():
    fig1 = plt.figure(1)
    #Plot Data-model
    frame1=fig1.add_axes((.1,.3,.8,.6))

def spectrum(x, y):
    plt.plot(x, y)
    
    plt.title("X-Ray Diffraction Pattern for Cu$_3$Au in Ordered Phase")
    plt.ylabel(count_axis)
    plt.xlabel(angle_axis)
    plt.show()

def peak(x_data, y_data, y_error, y_fit, bounds, popt, pcov, chi_sq, n):
    fig1 = plt.figure(1)
    #Plot Data-model
    frame1=fig1.add_axes((.1,.3,.8,.6))
    plt.errorbar(x_data, y_data, y_error, ls='none')
    plt.plot(x_data, y_fit)
    plt.ylabel(count_axis)
    plt.axis([bounds[0], bounds[1], -4, (abs(popt[3]) + abs(popt[1]))*1.5])
    
    plt.title("Peak Fit for Index " + str(n) + " in Cu$_3$Au Ordered Phase")
    
    chi_txt = str(int(chi_sq))
    ndf_txt = str(len(x_data) - 5)
    
    textstr = '$\mu='+value(abs(popt[0]), np.sqrt(pcov[0][0]))+'^\circ$\n' \
              '$A='+value(abs(popt[1]), np.sqrt(pcov[1][1]))+'$\n' \
              '$\sigma='+value(abs(popt[2]), np.sqrt(pcov[2][2]))+'^\circ$\n' \
              '$B='+value(abs(popt[3]), np.sqrt(pcov[3][3]))+'$\n' \
              '$\gamma='+value(abs(popt[4]), np.sqrt(pcov[4][4]))+'^\circ$\n' \
              '$\chi^2/\mathrm{ndf}='+chi_txt+' / '+ndf_txt+'$' 
    frame1.text(0.05, 0.95, textstr, transform=frame1.transAxes, fontsize=14,
        verticalalignment='top')
    
    frame1.set_xticklabels([]) #Remove x-tic labels for the first frame
    #plt.grid()
    
    frame2=fig1.add_axes((.1,.1,.8,.2))
    residual = (y_data - y_fit) / y_error
    res_error = y_error / y_error
    plt.errorbar(x_data, residual, res_error, ls='none')
    plt.plot(x_data, np.array([0]*len(x_data)))
    plt.axis([bounds[0], bounds[1], -4, 4 - 0.1])
    plt.xlabel(angle_axis)
    #plt.grid()
    
    plt.show()


def lattice(y_data, y_error, y_fit, x_data, popt, pcov, chi_sq):
    fig1 = plt.figure(1)
    #Plot Data-model
    frame1=fig1.add_axes((.1,.3,.8,.6))
    
    plt.errorbar(x_data, y_data*10, yerr=y_error*10, ls='none', fmt='-o')
    plt.plot(x_data, y_fit*10)
    plt.ylabel(lattice_axis)
    
    plt.title("Lattice Parameter Estimation for Cu$_3$Au Peaks in Ordered Phase")
    
    chi_txt = str(int(chi_sq))
    ndf_txt = str(len(x_data) - 2)
    
    textstr = '$a='+value(abs(popt[0])*10, np.sqrt(pcov[0][0])*10)+'\AA$\n' \
              '$\epsilon='+value(abs(popt[1]), np.sqrt(pcov[1][1]))+'$\n' \
              '$\chi^2/\mathrm{ndf}='+chi_txt+' / '+ndf_txt+'$'
    
    frame1.text(0.05, 0.95, textstr, transform=frame1.transAxes, fontsize=14,
        verticalalignment='top')
    
    frame1.set_xticklabels([]) #Remove x-tic labels for the first frame
    #plt.grid()
    
    frame2=fig1.add_axes((.1,.1,.8,.2))
    residual = (y_data - y_fit) / y_error
    res_error = y_error / y_error
    plt.errorbar(x_data, residual, res_error, ls='none', fmt='-o')
    plt.plot(x_data, np.array([0]*len(x_data)))
    plt.xlabel("2$\Theta$ ($^\circ$)")
    #plt.grid()
    
    plt.show()

def vegard(y, y_err, x, ni, cu):
    x_fit = np.array([x[0], x[-1]])
    y_fit = np.array([ni*10, cu*10])
    
    plt.errorbar(x, y*10, y_err*10)
    plt.plot(x_fit, y_fit)
    
    plt.ylabel(lattice_axis)
    plt.xlabel(composition_axis)
    plt.title("Lattice Parameter for Copper-Nickel Alloys")
    plt.show()

