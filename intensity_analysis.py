from __future__ import division
import parameter_analysis as p
import data as d
import math
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt


def interpolate(points, value):
    a = points[0]
    b = points[1]
    for point, i in zip(points, range(len(points) - 1)):
        if points[i + 1][0] > value:
            a = point
            b = points[i + 1]
    
    return a[1] + (b[1] - a[1])* (value - a[0])/(b[0] - a[0])

class IntensityAnalysis(p.ParameterAnalysis):
    def __init__(self, file_name=None, peaks=None, noise=None, uxd=False):
        p.ParameterAnalysis.__init__(self, file_name, peaks, noise, uxd)
        
        self.peak_ignore = 6
        
        self.peak_intensities = []
    
    def fill_intensities(self):
        for i in range(len(self.peaks)):
            self.peak_intensities.append(self.integrated_intensity(i))
    
    def cu_scattering_factor(self, angle):
        value = math.sin(d.rad(angle/2))/(self.x_ray_wavelength * 10)
        points = [(0.0, 29), (0.1, 25.9), (0.2, 21.6), (0.3, 17.9), (0.4, 15.2), (0.5, 13.3), (0.6, 11.7), (0.7, 10.2)]
        
        return interpolate(points, value)
        
    
    def au_scattering_factor(self, angle):
        value = math.sin(d.rad(angle/2))/(self.x_ray_wavelength * 10)
        points = [(0.0, 79), (0.1, 73.6), (0.2, 65.0), (0.3, 57.0), (0.4, 49.7), (0.5, 43.8), (0.6, 39.8), (0.7, 36.2)]
        
        return interpolate(points, value)
    
    def structure_factor(self, indices, order_parameter, angle):
        if p.all_even(indices) or p.all_odd(indices):
            return self.au_scattering_factor(angle) + 3 * self.cu_scattering_factor(angle)
        return order_parameter * (self.au_scattering_factor(angle) - self.cu_scattering_factor(angle))
        
    
    def intensity(self, i, order_parameter=1.0, epsilon=0.0):
        indices, _, __ = self.miller_indices("any")[i]
        p = [6, 12, 8, 6, 24, 24, 12, 30, 24, 24, 8, 24, 48][i]
        angle, error = self.peak_centres[i]
        #angle = angle - epsilon*2*180/np.pi*np.cos(d.rad(angle/2))
        TF = np.exp(epsilon * (np.sin(d.rad(angle/2)))**2)
        #if i == 6:
        #    angle = angle - 2
        lorentz_polarization = (1 + math.cos(d.rad(angle))**2) / (math.sin(d.rad(angle/2))**2 * math.cos(d.rad(angle/2)))
        
        return self.structure_factor(indices, order_parameter, angle)**2 * p * lorentz_polarization * TF
    
    def order_fit(self, x, order_parameter, epsilon):
        intensities = np.array([self.intensity(i + self.peak_ignore, order_parameter, epsilon) for i in range(len(x))])
        
        return intensities / np.sum(intensities)
    
    def compare_intensity(self, draw=False):
        x_data = []
        y_data = []
        y_error = []
        
        for i in range(len(self.peak_intensities)):
            x_data.append(self.peak_centres[i][0])
            y_data.append(self.peak_intensities[i][0])
            y_error.append(self.peak_intensities[i][1])
        
        x = np.array(x_data[self.peak_ignore:])
        y = np.array(y_data[self.peak_ignore:])
        y_err = np.array(y_error[self.peak_ignore:])
        
        y_err = y_err / np.sum(y)
        y = y / np.sum(y)
        
        popt, pcov = curve_fit(
            self.order_fit,
            x, y,
            sigma=y_err,
            p0=(1.0, 0.0)
        )
        
        if draw:
            print popt[0], popt[1]
            y_fit = self.order_fit(x, *popt)
            plt.errorbar(x, y, y_err)
            plt.plot(x, y_fit)
            plt.show()


if __name__ == "__main__":
    ###############################################
    peaks_A = [
        (500, 41.8, [41.5, 42.1]),
        (200, 48.6, [48.2, 49]),
        (140, 71.2, [70.8, 71.6]),
        (150, 86, [85.5, 86.5]),
        (35, 90.76, [90.3, 91.4])]
    noise_A = [87, 87.5]
    
    A = IntensityAnalysis(file_name="A_Peaks_3.txt", peaks=peaks_A, noise=noise_A)
    A.import_txt("MysteryA_LastPeak_2Feb.UXD", uxd=True)
    A.remove_noise()
    A.fill_peaks()
    
    A.average_lattice()
    
    #A.fill_peaks()
    #A.fill_intensities()
    
    #A.compare_intensity(draw=True)
    #A.average_lattice(draw=True)
    
    #A.draw()
    #for i in range(5):
    #    A.fit_double(i, draw=True, width=0.1)
    
    
    ###############################################
    peaks_B = [
        (120, 23.8, [23.3, 24.3]),
        (140, 33.9, [33.5, 34.3]),
        (600, 41.8, [41.4, 42.2]),
        (200, 48.7, [48.2, 49.2]),
        (60, 54.8, [54.2, 55.5]),
        (45, 60.6, [60, 61.2]),
        (60, 71.3, [70.8, 71.8]),
        (25, 76.4, [75.5, 77]),
        (15, 81.3, [80.5, 82]),
        (140, 86.2, [85.5, 86.7]),
        (40, 91, [90.3, 91.7]),
        (12, 95.9, [95, 96.5]),
        (20, 100.8, [99.5, 101.5])]
    
    noise_B = ([101.5, 102], [22, 22.5])
    
    B = IntensityAnalysis(file_name="MysteryB_31Jan.UXD", peaks=peaks_B, noise=noise_B, uxd=True)
    B.remove_noise()
    #d.Data("Old_data/Mystery_B_13-1-17.UXD").draw()
    #B.draw()
    B.fill_peaks()#draw_peaks=[1, 10])
    B.fill_intensities()
    
    B.average_lattice(draw=False, miller="any")
    B.compare_intensity(draw=True)
    
    #B.average_lattice(draw=True, miller="any")
    
    #B.draw()
    #for i in range(0, 1):
    #    B.fit_double(i, draw=True, width=0.1)
    
    #B.fill_peaks()
    
