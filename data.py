from __future__ import division
import numpy as np
import math
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

#TODO: Add second fit to funciton rather than lattice spacing / better sigma estimation from width

x_ray_wavelength = 0.1541

def mask(array, boolean_mask):
    temp_array = []
    for element, boolean in zip(array, boolean_mask):
        if boolean:
            temp_array.append(element)
    
    return np.array(temp_array)

def gaussian(x, a, b, c):
    return a * np.exp(-c * (x - b)*(x - b))

def linear(x, a, b):
    return a * x + b

def correction(x, a, b):
    return a * (1 - np.cos(x)**2 / np.sin(x) * b)

def chi_sq(y, y_fit, err):
    return sum((y - y_fit)**2 / err**2)/len(y)

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
        self.x_ray_wavelength = 0.1541

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
            p0=(peak, centre, 1.0/(width)**2)
        )
        
        y_fit = gaussian(x_data, *popt)
        chi_square = chi_sq(y_data, y_fit, y_error)
        print chi_square
        
        return popt, pcov
    
    def draw_peak(self, i):
        (peak, centre, width) = self.peaks[i]
        x_data, y_data, y_error = self.peak_data(i)
        params, _ = self.fit_gaussian(i)
        
        y_fit = gaussian(x_data, *params)
        
        plt.plot(x_data, y_data)
        plt.plot(x_data, y_fit)
        plt.axis([centre - width, centre + width, 0, peak*1.3])
        plt.show()
    
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
        error = abs(self.x_ray_wavelength*mode/2/math.sin(theta)**2*math.cos(theta)*theta_err)
        
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
        chi_square = chi_sq(y, y_fit, y_err)
        print chi_square
        
        if draw:
            plt.errorbar(thetas, y, yerr=y_err)
            plt.plot(thetas, y_fit)
            plt.show()
        
        print popt[0], np.sqrt(pcov[0][0]),  "\n"
        
        return popt[0], np.sqrt(pcov[0][0])
