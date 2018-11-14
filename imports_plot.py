# Common imports for plot

# Adjust plot defaults.
import matplotlib as mpl
# see mpl.rcParams for full list
mpl.rc('figure', figsize=(1600/240.0, 1200/240.0), dpi=240)
mpl.rc('lines', markersize=2.0)
mpl.rcParams['axes.formatter.limits'] = [-3,3]

import os
if "DISPLAY" not in os.environ:
    # When no X server, use this backend.
    mpl.use('Agg')

# For drawing
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure, clf, plot, scatter, xlabel, ylabel, xlim, ylim, matshow, colorbar, title, legend, semilogy, semilogx

plt.ion()

