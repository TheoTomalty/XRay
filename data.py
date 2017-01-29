from __future__ import division
import numpy as np
import math
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

#TODO: Add second fit to funciton rather than lattice spacing / better sigma estimation from width

x_ray_wavelength = 0.15418
x_ray_separation  = 0.0004

statistical_xray = x_ray_separation/x_ray_wavelength/20

def mask(array, boolean_mask):
    temp_array = []
    for element, boolean in zip(array, boolean_mask):
        if boolean:
            temp_array.append(element)
    
    return np.array(temp_array)

def gaussian(x, a, b, c):
    return a * np.exp(-(x - b)*(x - b)/2/c**2)

def linear(x, a, b):
    return a * x + b

def correction(x, a, b):
    return a * (1 - np.cos(x)**2 / np.sin(x) * b)

def chi_sq(y, y_fit, err, num_params):
    return sum((y - y_fit)**2 / err**2)/(len(y) - num_params)

def all_even(elements):
    for element in elements:
        if not element % 2:
            return False
    
    return True

def all_odd(elements):
    for element in elements:
        if element % 2:
            return False
    
    return True

class Data(object):
    def __init__(self, file_name=None):
        self.count = np.array([])
        self.angle = np.array([])
        self.error = np.array([])
        
        self.peaks = None
        
        self.noise_cut = 40
        self.angle_step = 0.05
        self.x_ray_wavelength = x_ray_wavelength

        if file_name is not None:
            self.import_txt(file_name)

    def import_txt(self, file_name):
        with open(file_name, "r") as file:
            for line in file:
                angle, count = map(float, line[1:].split(","))
                self.count = np.append(self.count, count)
                self.error = np.append(self.error, np.sqrt(count))
                self.angle = np.append(self.angle, angle)
        
        self.peaks = self.peak_finder()
    
    def draw(self):
        plt.plot(self.angle, self.count)
        plt.title('X-Ray Diffraction Spectrum of Copper')
        plt.ylabel('Photon Count')
        plt.xlabel('2$\Theta$ ($^\circ$)')
        plt.show()
    
    def index(self, angle):
        index_range = len(self.angle)
        angle_begin = self.angle[0]
        angle_end = self.angle[-1]
        
        ratio = (angle - angle_begin)/(angle_end - angle_begin)
        
        return round(ratio * index_range)
    
    def peak_finder(self):
        parameters = []
        
        start = None
        step_y = None
        step_x = None
        
        for count, angle in zip(self.count, self.angle):
            if count > self.noise_cut:
                if step_y is None:
                    start = angle
                    step_y = count
                    step_x = angle
                elif count > step_y:
                    step_y = count
                    step_x = angle
            elif start is not None:
                width = angle - start
                if width > self.angle_step*7:
                    parameters.append((step_y, step_x, width*2))
                
                start = None
                step_y = None
                step_x = None
        
        return parameters
    
    def peak_data(self, i):
        (peak, centre, width) = self.peaks[i]
        data_range = [self.index(centre - width/2), self.index(centre + width/2)]
        
        x_data = self.angle[data_range[0]:data_range[1]]
        y_data = self.count[data_range[0]:data_range[1]]
        y_error = self.error[data_range[0]:data_range[1]]
        
        boolean = [y > peak/4 for y in y_data]
        
        return mask(x_data, boolean), mask(y_data, boolean), mask(y_error, boolean)
    
    def fit_gaussian(self, i):
        (peak, centre, width) = self.peaks[i]
        x_data, y_data, y_error = self.peak_data(i)
        
        popt, pcov = curve_fit(
            gaussian,
             x_data, y_data,
            sigma=y_error,  
            p0=(peak, centre, width/2)
        )
        
        y_fit = gaussian(x_data, *popt)
        chi_square = chi_sq(y_data, y_fit, y_error, 3)
        print "Gaussian", popt[1], np.sqrt(pcov[1][1]), chi_square*(len(y_data) - 3), len(y_data) - 3
        
        return popt, pcov
    
    def draw_peak(self, i):
        (peak, centre, width) = self.peaks[i]
        x_data, y_data, y_error = self.peak_data(i)
        params, _ = self.fit_gaussian(i)
        
        y_fit = gaussian(x_data, *params)
        
        plt.errorbar(x_data, y_data, y_error, ls='none', fmt='-o')
        
        x_cont = np.arange(x_data[0], x_data[-1], 0.005)
        y_cont = gaussian(x_cont, *params)
        plt.plot(x_cont, y_cont)
        plt.ylabel("Photon Count")
        width1 = x_data[-1] - x_data[0]
        range = [centre - width1, centre + width1, 0, peak*1.2]
        plt.axis(range)
        
        #print "Fit Params: ", params, np.sqrt(_)
        
        return y_data, y_error, y_fit, x_data, range
    
    def milner_indices(self):
        milner = []
        test_max = 5
        for i in range(test_max):
            for j in range(test_max):
                for k in range(test_max):
                    n = (i, j, k)
                    if all_even(n) or all_odd(n):
                        sq = n[0]**2 + n[1]**2 + n[2]**2
                        milner.append((n, sq))
        
        milner = sorted(milner, key=lambda key: key[1])
        index = 0
        del milner[0]
        while index < len(milner) - 1:
            if milner[index][1] == milner[index + 1][1]:
                del milner[index]
            else:
                index += 1
        
        return milner
    
    def peak_centre(self, i):
        popt, pcov = self.fit_gaussian(i)
        
        return popt[1], np.sqrt(pcov[1][1])
    
    def lattice_parameter(self, i):
        centre, cent_err = self.peak_centre(i)
        
        theta = math.radians(centre/2)
        theta_err = math.radians(cent_err/2)
        mode = np.sqrt(self.milner_indices()[i][1])
        
        spacing = self.x_ray_wavelength*mode/(2*math.sin(theta))
        error1 = abs(self.x_ray_wavelength*mode/2/math.sin(theta)**2*math.cos(theta)*theta_err)
        error2 = statistical_xray * spacing
        error = math.sqrt(error1**2 + error2**2)
        
        print "Lattice Parameter", spacing, error
        
        return spacing, error, centre
    
    def average_lattice(self, draw=False):
        pts = [self.lattice_parameter(i) for i in range(len(self.peaks))]
        y = np.array([pt[0] for pt in pts])
        y_err = np.array([pt[1] for pt in pts])
        x = [pt[2] for pt in pts]
        
        thetas = np.array([math.radians(point)/2 for point in x])
        
        popt, pcov = curve_fit(
            correction,
            thetas, y,
            sigma=y_err,
            p0=(0.1/100.0, 0.36)
        )
        
        y_fit = correction(thetas, *popt)
        chi_square = chi_sq(y, y_fit, y_err, 2)
        #print chi_square*(len(y) - 2), (len(y) - 2)
        
        if draw:
            plt.errorbar(x, y*10, yerr=y_err*10, ls='none', fmt='-o')
            x_cont = np.arange(x[0]-1, x[-1] + 2, 1)
            theta_cont = np.array([math.radians(point)/2 for point in x_cont])
            y_cont = correction(theta_cont, *popt)
            plt.plot(x_cont, y_cont*10)
            plt.ylabel("Lattice Parameter ($\AA$)")
        
        #print "Combined Lattice", popt[0], np.sqrt(pcov[0][0]), popt[0]*x_ray_separation/x_ray_wavelength/2, "\n"
        #print popt, np.sqrt(pcov)
        
        return y, y_err, y_fit, x
