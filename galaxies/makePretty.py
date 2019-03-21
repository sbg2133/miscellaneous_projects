import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
plt.ion()
plt.rcParams['axes.grid'] = True

def pretty():
    ax = plt.gca()
    ax.xaxis.set_major_locator(MaxNLocator(prune='lower'))
    [i.set_linewidth(3.) for i in ax.spines.itervalues()]
    ax.tick_params(axis='x', labelsize=18)
    ax.tick_params(axis='y', labelsize=18)
    [i.set_linewidth(2.) for i in ax.spines.itervalues()]
    return
