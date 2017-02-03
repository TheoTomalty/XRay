from __future__ import division
import numpy as np
import math
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

#TODO: Add second fit to funciton rather than lattice spacing / better sigma estimation from width

k_alpha1 = 0.154056
k_alpha2 = 0.154439
x_ray_wavelength = 0.15418
x_ray_separation  = k_alpha2 - k_alpha1

a = (x_ray_wavelength - k_alpha2)/(-x_ray_separation)
b = (x_ray_wavelength - k_alpha1)/(x_ray_separation)

proportion = b/a

statistical_xray = 0

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

#def double_peak(x, centre, A, gamma1, gamma2, B):
#    centre1 = centre + 2*deg((k_alpha1/x_ray_wavelength - 1) * np.tan(rad(centre/2)))
#    centre2 = centre + 2*deg((k_alpha2/x_ray_wavelength - 1) * np.tan(rad(centre/2)))
#    
#    scaling1 = (x - centre1)/gamma1
#    scaling2 = (x - centre2)/gamma2
#    
#    return A * (1/(1 + scaling1**2) + (proportion * gamma1/gamma2) / (1 + scaling2**2)) + B

def double_peak(x, centre, A, sigma1, sigma2, tail, gamma):
    centre1 = centre + 2*deg((k_alpha1/x_ray_wavelength - 1) * np.tan(rad(centre/2)))
    centre2 = centre + 2*deg((k_alpha2/x_ray_wavelength - 1) * np.tan(rad(centre/2)))
    
    scaling1 = -(x - centre1)**2 / (2*sigma1**2)
    scaling2 = -(x - centre2)**2 / (2*sigma1**2)
    lorentz1 = (x - centre1)**2/gamma**2
    lorentz2 = (x - centre2)**2/gamma**2
    
    return A * (np.exp(scaling1) + proportion * np.exp(scaling2))\
           + A*tail * (1/(1 + lorentz1) + proportion / (1 + lorentz2))

def noise_funct(x, A):
    return np.array([A]*len(x)) #+ B * np.exp(-x**2/2/sigma**2)

def linear(x, a, b):
    return a * x + b

def chi_sq(y, y_fit, err, num_params):
    return sum((y - y_fit)**2 / err**2)/(len(y) - num_params)


class Data(object):
    def __init__(self, file_name=None, peaks=None, noise=None, uxd=False):
        self.count = np.array([])
        self.angle = np.array([])
        self.error = np.array([])
        
        self.peaks = peaks
        self.noise = noise
        
        self.noise_cut = 30
        self.angle_step = 0.05
        self.x_ray_wavelength = x_ray_wavelength

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
                self.error = np.append(self.error, np.sqrt(count))
                self.angle = np.append(self.angle, angle)
    
    def draw(self):
        plt.plot(self.angle, self.count)
        #plt.show()
    
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
    
    def remove_noise(self, draw=False):
        if self.noise == None:
            boolean = [True]*len(self.angle)
            ranges = []
            
            for peak, centre, width in self.peaks:
                ranges.append((self.index(centre - width), self.index(centre + width)))
            
            for bounds in ranges:
                for index in range(bounds[0], bounds[1] + 1):
                    boolean[index] = False
            
            x = mask(self.angle, boolean)
            y = mask(self.count, boolean)
        else:
            x, y, _ = self.get_range(self.noise[0], self.noise[1])
        
        average = np.mean(y)
        rms = np.sqrt(np.mean((y - np.array([average]*len(x)))**2))
        
        if draw:
            x_cont = np.arange(x[0], x[-1], 0.005)
            y_cont = np.array([average]*len(x_cont))
            
            plt.plot(x, y)
            plt.plot(x_cont, y_cont)
            plt.show()
        
        count_correction = np.array([average]*len(self.angle))
        error_correction = np.array([rms]*len(self.angle))
        
        self.count -= count_correction
        self.error = np.sqrt(error_correction**2 + self.error**2)
    
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
    
    def fit_double(self, i, draw=False):
        peak, centre, bounds = self.peaks[i]
        x_data, y_data, y_error = self.get_range(bounds[0], bounds[1])
        
        popt, pcov = curve_fit(
            double_peak,
             x_data, y_data,
            sigma=y_error,  
            p0=(centre, peak, 0.1, 0.1, 0, (bounds[0] - bounds[1])/4)
        )
        
        if draw:
            print popt[4]
            plt.plot(x_data, y_data)
            y_fit = double_peak(x_data, *popt)
            plt.plot(x_data, y_fit)
            plt.show()
        
        return popt, pcov
    
    def fit_gaussian(self, i, draw=False):
        peak, centre, width = self.peaks[i]
        x_data, y_data, y_error = self.get_range(centre - width/2, centre + width/2)
        if len(x_data) < 10 or peak < 100:
            print "ERROR", len(x_data), peak
            return None, None
        
        popt, pcov = curve_fit(
            gaussian,
             x_data, y_data,
            sigma=y_error,  
            p0=(peak/4, centre, width/4, peak, width/6)
        )
        
        y_fit = gaussian(x_data, *popt)
        chi_square = chi_sq(y_data, y_fit, y_error, 5)
        print "Gaussian", popt[1], np.sqrt(pcov[1][1]), chi_square*(len(y_data) - 3), len(y_data) - 3
        
        if draw:
            plt.errorbar(x_data, y_data, y_error, ls='none', fmt='-o')
            
            x_cont = np.arange(x_data[0], x_data[-1], 0.005)
            y_cont = gaussian(x_cont, *popt)
            plt.plot(x_cont, y_cont)
            plt.show()
        
        return popt, pcov
    
    def integrated_intensity(self, i, draw=False):
        peak, centre, width = self.peaks[i]
        x_data, y_data, y_error = self.get_range(centre - width/2, centre + width/2)
        if len(x_data) < 10 or peak < 100:
            print "ERROR", len(x_data), peak
            return None, None
            
        integral = 0
        variance = 0
        for val, error in zip(y_data, y_error):
            integral += val * self.angle_step
            variance += (error * self.angle_step)**2
        
        if draw:
            plt.errorbar(x_data, y_data, y_error)
            plt.show()
        
        return integral, np.sqrt(variance)
