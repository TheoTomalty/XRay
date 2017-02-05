import data as d
import parameter_analysis as p
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

#TODO: Make B peaks better fits by removing tails

class MultiLattice(d.Data):
    def __init__(self, file_name=None, peaks=None, noise=None, uxd=False, draw=False, initialize_sn=False, initialize_pb=False):
        d.Data.__init__(self, file_name, None, noise, uxd)
        self.miller_sn = [
            (2, 0, 0), (1, 0, 1), (2, 2, 0), (2, 1, 1), (3, 0, 1), 
            (1, 1, 2), (4, 0, 0), (3, 2, 1), (4, 2, 0), (4, 1, 1), 
            (3, 1, 2), (4, 3, 1)
        ]
        self.global_factor = []
        
        pb_peaks = []
        sn_peaks = []
        sn_counter = 0
        for peak, centre, bounds, type in peaks:
            if type == "Sn":
                sn_counter += 1
                sn_peaks.append((peak, centre, bounds))
            elif type == "Pb":
                pb_peaks.append((peak, centre, bounds))
            elif type == "Sn_N":
                del self.miller_sn[sn_counter]
            
        self.Sn = p.ParameterAnalysis(file_name, sn_peaks, noise, uxd)
        self.Pb = p.ParameterAnalysis(file_name, pb_peaks, noise, uxd)
        
        if initialize_sn:
            if noise is not None:
                self.Sn.remove_noise()
            self.Sn.fill_peaks(draw=draw)
            self.a, self.epsilon = self.lattice_a(draw=draw)
        if initialize_pb:
            if noise is not None:
                self.Pb.remove_noise()
            self.Pb.fill_peaks(draw=draw)
    
    def lattice_a(self, draw=False):
        points = []
        for indices, i in zip(self.miller_sn, range(len(self.miller_sn))):
            if indices[2] == 0:
                x, dx = self.Sn.peak_centres[i]
                theta, d_theta = d.rad(x/2), d.rad(dx/2)
                param = self.Sn.x_ray_wavelength/(2 * np.sin(theta))
                y = np.sqrt(indices[0]**2 + indices[1]**2) * param
                y_err = abs(y * np.cos(theta) / np.sin(theta) * d_theta)
                points.append((y, y_err, x))
        
        print points
        return self.Sn.average_lattice(pts=points, draw=draw)
    
    def correction(self, x, b, epsilon):
        factors = np.array(self.global_factor)
        return b * (1 - factors * np.cos(x) / np.sin(x) * epsilon)
    
    def lattice_b(self, draw=False, remove_peaks=None):
        pts = []
        sum = 0
        variance = 0
        for indices, i in zip(self.miller_sn, range(len(self.miller_sn))):
            if indices[2] != 0:
                x, dx = self.Sn.peak_centres[i]
                theta, d_theta = d.rad(x/2), d.rad(dx/2)
                param = self.Sn.x_ray_wavelength/(2 * np.sin(theta))
                C = (indices[0]**2 + indices[1]**2)/self.a**2
                
                y = indices[2]/np.sqrt(1/param**2 - C)
                y_err = abs(y/(1 - param**2 * C) * np.cos(theta) / np.sin(theta) * d_theta)
                
                sum += y
                variance += y_err**2
                self.global_factor.append(1/(1 - param**2 * C))
                pts.append((y, y_err, x))
        
        if remove_peaks is not None:
            for num in remove_peaks:
                del pts[num]
                del self.global_factor[num]
        
        y1 = np.array([pt[0] for pt in pts])
        y_error = np.array([pt[1] for pt in pts])
        #y_error = np.array([0.0001]*len(y1))
        x1 = [pt[2] for pt in pts]
        
        thetas = np.array([d.rad(point/2) for point in x1])
        
        popt, pcov = curve_fit(
            self.correction,
            thetas, y1,
            sigma=y_error,
            p0=(0.31750, 0.0004)
        )
        print popt
        if draw:
            y_fit = self.correction(thetas, *popt)
            plt.errorbar(x1, y1*10, yerr=y_error*10, ls='none', fmt='-o')
            plt.plot(x1, y_fit*10)
            plt.show()

peaks_Pb0 = [
    (250, 30.7, [30.4, 31], "Sn"),
    (350, 32.1, [31.8, 32.4], "Sn"),
    (140, 43.95, [43.6, 44.4], "Sn"),
    (250, 45, [44.6, 45.4], "Sn"),
    (50, 55.4, [55, 56], "Sn"),
    (90, 62.6, [62.2, 63.2], "Sn"),
    (15, 63.9, [63.4, 64.3], "Sn"),
    (60, 64.65, [64.4, 65.1], "Sn"),
    (65, 72.5, [72.2, 73], "Sn"),
    (40, 73.2, [73, 73.7], "Sn"),
    (65, 79.6, [79.2, 80.1], "Sn"),
    (150, 89.45, [89.2, 90], "Sn")
]

noise_Pb0=([32.7, 33], [90.1, 90.5])

Pb0 = MultiLattice("Pb0_2Feb.UXD", peaks=peaks_Pb0, noise=noise_Pb0, uxd=True, initialize_sn=True)

#for i in range(len(peaks)):
#    Pb0.fit_double(i, width=0.1, draw=True)

peaks_Pb100 = [
    (1200, 31.4, [31, 31.8], "Pb"),
    (500, 36.4, [36, 36.8], "Pb"),
    (240, 52.4, [51.8, 52.8], "Pb"),
    (260, 62.2, [61.8, 62.8], "Pb"),
    (120, 65.3, [64.8, 65.8], "Pb"),
    (30, 77.2, [76.5, 78], "Pb"),
    (80, 85.5, [85.2, 86.1], "Pb"),
    (80, 88.3, [87.8, 89], "Pb")
]
noise_Pb100 = ([51, 51.5], [86.5, 87])

Pb100 = MultiLattice("Pb100_30_1_17.UXD", peaks=peaks_Pb100, noise=noise_Pb100, uxd=True, initialize_pb=True)
#Pb100.Pb.average_lattice(draw=True)
#Pb100.draw()

peaks_Pb25 = sorted(
    [(0.25*peak[0], peak[1], peak[2], peak[3]) for peak in peaks_Pb0]
    + [(0.75*peak[0], peak[1], peak[2], peak[3]) for peak in peaks_Pb100],
    key=lambda key: key[1]
)
peaks_Pb25[8] = (peaks_Pb25[8][0], 62.25, [61.8, 62.5], 'Pb')
peaks_Pb25[9] = (peaks_Pb25[9][0], 62.6, [62.5, 63.2], 'Sn')
peaks_Pb25[12] = (peaks_Pb25[12][0], 65.4, [65.1, 65.8], 'Pb')
#print peaks_middle

noise_Pb25 = ([29, 29.5], [98, 99])
Pb25 = MultiLattice("Pb25_1Feb.UXD", peaks=peaks_Pb25, noise=noise_Pb25, uxd=True, initialize_sn=True, initialize_pb=True)
#Pb25.lattice_b(draw=True, remove_peaks=[-1, -3])

###############################################################

peaks_Pb50 = sorted(
    [(0.5*peak[0], peak[1], peak[2], peak[3]) for peak in peaks_Pb0]
    + [(0.5*peak[0], peak[1], peak[2], peak[3]) for peak in peaks_Pb100],
    key=lambda key: key[1]
)
peaks_Pb50[8] = (peaks_Pb50[8][0], 62.25, [61.8, 63], 'Pb')
peaks_Pb50[9] = (0, 0, [0, 0], 'Sn_N')
peaks_Pb50[12] = (peaks_Pb50[12][0], 65.4, [65.1, 65.8], 'Pb')
#print peaks_Pb50

noise_Pb50 = [86.6, 87.6]
Pb50 = MultiLattice("Pb_50_02_01_17.UXD", peaks=peaks_Pb50, uxd=True, initialize_sn=True, initialize_pb=True)
#Pb50.draw()
#Pb50.lattice_b(draw=True, remove_peaks=[-1, -2, -2])


###############################################################

peaks_Pb75 = sorted(
    [(0.25*peak[0], peak[1], peak[2], peak[3]) for peak in peaks_Pb0]
    + [(0.75*peak[0], peak[1], peak[2], peak[3]) for peak in peaks_Pb100],
    key=lambda key: key[1]
)
#peaks_Pb75[2] = (40, 32.15, [31.9, 32.4], peaks_Pb75[2][3])
peaks_Pb75[4] = (0, 0, [0, 0], 'Sn_N')
peaks_Pb75[7] = (0, 0, [0, 0], 'Sn_N')
peaks_Pb75[8] = (peaks_Pb75[8][0], 62.25, [61.8, 63], 'Pb')
peaks_Pb75[10] = (0, 0, [0, 0], 'Sn_N')
peaks_Pb75[13] = (0, 0, [0, 0], 'Sn_N')
peaks_Pb75[14] = (0, 0, [0, 0], 'Sn_N')
peaks_Pb75[16] = (0, 0, [0, 0], 'Sn_N')
peaks_Pb75[-1] = (0, 0, [0, 0], 'Sn_N')
peaks_Pb75[12] = (peaks_Pb75[12][0], 65.4, [65.1, 65.8], 'Pb')
print peaks_Pb75

noise_Pb75 = ([53, 53.5], [107, 107.5])
Pb75 = MultiLattice("Pb75_31Jan.UXD", peaks=peaks_Pb75, noise=noise_Pb75, uxd=True, initialize_pb=True)

Pb0.lattice_b(draw=True, remove_peaks=[-1, -2])
Pb25.lattice_b(remove_peaks=[-1, -3])
Pb50.lattice_b(remove_peaks=[-1, -2, -2])

Pb25.Pb.average_lattice(draw=True)
Pb50.Pb.average_lattice(draw=True)
Pb75.Pb.average_lattice(draw=True)
Pb100.Pb.average_lattice(draw=True)
