import data as d
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

###############################
peaks = [
    (500, 41.8, [40.5, 43]),
    (200, 48.6, [48, 49.5]),
    (140, 71.2, [70.5, 72]),
    (150, 86, [85, 87]),
    (30, 90.9,[90.3, 91.4])]
noise = [64, 67]

A = d.Data(file_name="A_Peaks_3.txt", peaks=peaks, noise=noise)
A.import_txt("MysteryA_LastPeak_2Feb.UXD", uxd=True)
A.remove_noise()

#for i in range(5):
#    A.fit_double(i, draw=True)

###############################
peaksB = [
    (120, 23.8, [23, 24.5]),
    (140, 33.9, [33, 34.5]),
    (600, 41.8, [41, 42.5]),
    (200, 48.7, [47.5, 49.5]),
    (60, 54.8, [54, 56]),
    (45, 60.6, [59.5, 61.5]),
    (60, 71.3, [70.5, 72]),
    (25, 76.4, [75.5, 77.5]),
    (15, 81.3, [80.5, 82]),
    (140, 86.2, [85, 87]),
    (40, 91, [90, 92]),
    (12, 95.9, [95, 96.5]),
    (20, 100.8, [99.5, 102])]
noiseB = [22, 23]

B = d.Data(file_name="MysteryB_31Jan.UXD", peaks=peaksB, noise=noiseB, uxd=True)
B.remove_noise()

#for i in range(13):
#    B.fit_double(i, draw=True)
