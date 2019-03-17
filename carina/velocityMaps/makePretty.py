import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
plt.ion()
#plt.rcParams['axes.grid'] = True
matplotlib.rcParams.update({'font.size': 18})

def pretty():
    ax = plt.gca()
    ax.xaxis.set_major_locator(MaxNLocator(prune='lower'))
    [i.set_linewidth(2.) for i in ax.spines.itervalues()]
    ax.tick_params(axis='x', labelsize=18,\
          bottom = True, top = True, left = True, right = True,\
           width = 2, direction = 'in')
    ax.tick_params(axis='y', labelsize=15,\
          bottom = True, top = True, left = True, right = True,\
           width = 2, direction = 'in')
    [i.set_linewidth(2.) for i in ax.spines.itervalues()]
    return
