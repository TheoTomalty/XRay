import data as d
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

#d.Data("Pb0_2Feb.UXD", uxd=True).draw()
#d.Data("Pb25_1Feb.UXD", uxd=True).draw()
#d.Data("Pb_50_02_01_17.UXD", uxd=True).draw()
#d.Data("Pb75_31Jan.UXD", uxd=True).draw()
#d.Data("Pb100_30_1_17.UXD", uxd=True).draw()
#
#plt.show()

D1 = d.Data("Cu_0_Peaks.UXD", uxd=True).draw()
#D2 = d.Data("Cu_25_26_1_17.UXD", uxd=True).draw()
#D3 = d.Data("Cu_50_Peaks.UXD", uxd=True).draw()
#D4 = d.Data("Cu_75_Peaks.UXD", uxd=True).draw()
#
#D5 = d.Data("Cu_Peak1.UXD", uxd=True)
#D5.import_txt("Cu_Peak2.UXD", uxd=True)
#D5.import_txt("Cu_Peak3.UXD", uxd=True)
#D5.import_txt("Cu_Peak4.UXD", uxd=True)
#D5.import_txt("Cu_Peak5.UXD", uxd=True)
#D5.draw()


plt.show()
