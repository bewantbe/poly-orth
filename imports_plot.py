# Common imports for plot

# Drawing
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure, clf, plot, scatter, xlabel, ylabel, xlim, ylim, matshow, colorbar, title, legend

# Adjust plot defaults.
import matplotlib as mpl
# see mpl.rcParams for full list
mpl.rc('figure', figsize=(1600/240.0, 1200/240.0), dpi=240)
mpl.rc('lines', markersize=2.0)
mpl.rcParams['axes.formatter.limits'] = [-3,3]
plt.ion()

