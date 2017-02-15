from __future__ import division
import numpy as np
import math
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import draw as display

k_alpha1 = 0.154056
k_alpha2 = 0.154439
x_ray_wavelength = 0.15418
x_ray_separation  = k_alpha2 - k_alpha1

a = (x_ray_wavelength - k_alpha2)/(-x_ray_separation)
b = (x_ray_wavelength - k_alpha1)/(x_ray_separation)

proportion = b/a

def rad(angle):
    return angle*math.pi/180

def deg(angle):
    return angle * 180/math.pi

def mask(array, boolean_mask):
    temp_array = []
    for element, boolean in zip(array, boolean_mask):
        if boolean:
            temp_array.append(element)
    
    return np.array(temp_array)

def gaussian(x, A, x_0, sigma, B, gamma):
    return A * np.exp(-(x - x_0)**2/2/sigma**2) + B/(np.pi*gamma*(1.0 + (x - x_0)**2/gamma**2))

def double_peak(x, centre, A, sigma1, tail, gamma):
    centre1 = centre
    centre2 = centre + 2*deg((x_ray_separation/x_ray_wavelength) * np.tan(rad(centre/2)))
    
    scaling1 = -(x - centre1)**2 / (2*sigma1**2)
    scaling2 = -(x - centre2)**2 / (2*sigma1**2)
    lorentz1 = (x - centre1)**2/gamma**2
    lorentz2 = (x - centre2)**2/gamma**2
    
    return A * (np.exp(scaling1) + proportion * np.exp(scaling2)) \
           + tail * (1/(1 + lorentz1) + proportion / (1 + lorentz2))

def noise_funct(x, A):
    return np.array([A]*len(x)) #+ B * np.exp(-x**2/2/sigma**2)

def linear(x, a, b):
    return a * x + b

def chi_sq(y, y_fit, err, num_params):
    return sum((y - y_fit)**2 / err**2)

class Data(object):
    def __init__(self, file_name=None, peaks=None, noise=None, uxd=False):
        self.count = np.array([])
        self.angle = np.array([])
        self.error = np.array([])
        
        self.peaks = peaks
        self.noise = noise
        
        self.noise_cut = 30
        self.angle_step = 0.05
        self.x_ray_wavelength = k_alpha1

        if file_name is not None:
            self.import_txt(file_name, uxd)
            if self.peaks == None:
                self.peaks = self.peak_finder()
    
    def import_txt(self, file_name, uxd=False):
        if uxd:
            min = 1
            max = -1
            string = "      "
        else:
            min = 1
            max = -1
            string = ","
        with open(file_name, "r") as file:
            for line in file:
                angle, count = map(float, line[min:max].split(string))
                self.count = np.append(self.count, count)
                self.error = np.append(self.error, (np.sqrt(count) if count else 1))
                self.angle = np.append(self.angle, angle)
    
    def draw(self):
        display.spectrum(self.angle, self.count)
    
    def index(self, angle):
        index_range = len(self.angle)
        angle_begin = self.angle[0]
        angle_end = self.angle[-1]
        
        ratio = (angle - angle_begin)/(angle_end - angle_begin)
        
        return int(round(ratio * index_range))
    
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
                parameters.append((step_y, step_x, width*2))
                
                start = None
                step_y = None
                step_x = None
        
        return parameters
    
    def remove_noise(self):
        if isinstance(self.noise, tuple):
            begin_bounds, end_bounds = self.noise
            x1, y1, _ = self.get_range(begin_bounds[0], begin_bounds[1])
            x2, y2, _ = self.get_range(end_bounds[0], end_bounds[1])
            
            centre_1, average_1 = np.mean(x1), np.mean(y1)
            centre_2, average_2 = np.mean(x2), np.mean(y2)
            
            count_correction = np.array(
                [average_1 + (average_2 - average_1) * (angle - centre_1)/(centre_2 - centre_1) for angle in self.angle]
            )
            
        elif self.noise == None:
            boolean = [True]*len(self.angle)
            ranges = []
            
            for peak, centre, width in self.peaks:
                ranges.append((self.index(centre - width), self.index(centre + width)))
            
            for bounds in ranges:
                for index in range(bounds[0], bounds[1] + 1):
                    boolean[index] = False
            
            y = mask(self.count, boolean)
            count_correction = np.array([np.mean(y)]*len(self.angle))
            
        else:
            __, y, _ = self.get_range(self.noise[0], self.noise[1])
            count_correction = np.array([np.mean(y)]*len(self.angle))
        
        self.count -= count_correction
    
    def get_range(self, min, max):
        index_min = None
        index_max = 0
        for x, index in zip(self.angle, range(len(self.angle))):
            if index_min is None and x > min:
                index_min = index
            if x < max:
                index_max = index
        
        x_data = self.angle[index_min:index_max]
        y_data = self.count[index_min:index_max]
        y_error = self.error[index_min:index_max]
        
        return x_data, y_data, y_error
    
    def fit_double(self, i, draw=False, width=0.5):
        peak, centre, bounds = self.peaks[i]
        x_data, y_data, y_error = self.get_range(bounds[0], bounds[1])
        
        popt, pcov = curve_fit(
            double_peak,
             x_data, y_data,
            sigma=y_error,  
            p0=(centre, peak, width, 0, (bounds[1] - bounds[0])/2)
        )
        
        if draw:
            y_fit = double_peak(x_data, *popt)
            if i < 7:
                n = (1, 1, 0)
            else:
                n = (2, 2, 2)
            display.peak(x_data, y_data, y_error, y_fit, bounds, popt, pcov, chi_sq(y_data, y_fit, y_error, 5), n)
        
        x = popt[0]
        systematic_error = 1/15 * x_ray_separation/x_ray_wavelength * np.tan(rad(x/2)) * 180/np.pi
        
        return popt[0], np.sqrt(pcov[0][0]**2 + systematic_error**2)
    
    def integrated_intensity(self, i, draw=False):
        peak, centre, bounds = self.peaks[i]
        x_data, y_data, y_error = self.get_range(bounds[0] - 2, bounds[1] + 2)
        print (bounds[0] - 2, bounds[1] + 2)
        
        base = (x_data[-1] - x_data[0])
        
        integral = 0
        variance = 0
        for val, error in zip(y_data, y_error):
            integral += val * self.angle_step
            variance += (error * self.angle_step)**2
        
        if draw:
            plt.plot(x_data, y_data)
            plt.show()
        
        return integral, np.sqrt(variance + (3*base)**2)
