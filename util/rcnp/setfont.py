import matplotlib.pyplot as pylab

def setfont(font='helvetica',unicode=True):
    r"""
    Set Matplotlibs rcParams to use LaTeX for font rendering.
    Revert all changes by calling rcdefault() from matplotlib.

    Parameters:
    -----------
    font: string
        "Helvetica"
        "Times"
        "Computer Modern"

    usetex: Boolean
        Use unicode. Default: False.
    """

    # Use TeX for all figure text!
    pylab.rc('text', usetex=True)

    font = font.lower().replace(" ","")
    if font == 'times':
        # Times
        font = {'family':'serif', 'serif':['Times']}
        preamble  = r"""
                       \usepackage{color}
                       \usepackage{mathptmx}
                    """
    elif font == 'helvetica':
        # Helvetica
        # set serif, too. Otherwise setting to times and then
        # Helvetica causes an error.
        font = {'family':'sans-serif','sans-serif':['Helvetica'],
                'serif':['cm10']}
        preamble  = r"""
                       \usepackage{color}
                       \usepackage[tx]{sfmath}
                       \usepackage{helvet}
                    """
    else:
        # Computer modern serif
        font = {'family':'serif', 'serif':['cm10']}
        preamble  = r"""
                       \usepackage{color}
                    """

    if unicode:
        # Unicode for Tex
        #preamble =  r"""\usepackage[utf8]{inputenc}""" + preamble
        # inputenc should be set automatically
        pylab.rcParams['text.latex.unicode']=True

    #print font, preamble
    pylab.rc('font',**font)
    pylab.rcParams['text.latex.preamble'] = preamble


def setticks(ax,xlog=False,ylog=False,xmajor=5,xminor=1,ymajor=2,yminor=0.5):

    if not xlog:
        xmajorLocator   = pylab.MultipleLocator(xmajor)
        xmajorFormatter = pylab.FormatStrFormatter('%d')
        xminorLocator   = pylab.MultipleLocator(xminor)
        ax.xaxis.set_major_locator(xmajorLocator)
        #ax.xaxis.set_major_formatter(xmajorFormatter)
        ax.xaxis.set_minor_locator(xminorLocator)

    if not ylog:
        ymajorLocator   = pylab.MultipleLocator(ymajor)
        ymajorFormatter = pylab.FormatStrFormatter('%d')
        yminorLocator   = pylab.MultipleLocator(yminor)
        ax.yaxis.set_major_locator(ymajorLocator)
        #ax.yaxis.set_major_formatter(ymajorFormatter)
        ax.yaxis.set_minor_locator(yminorLocator)

    ax.get_yaxis().set_tick_params(which='both', direction='out')
    ax.get_xaxis().set_tick_params(which='both', direction='out')
    for tick in ax.xaxis.get_ticklines():
        tick.set_markersize(4.5)
    for tick in ax.yaxis.get_ticklines():
        tick.set_markersize(4.5)
    for tick in ax.xaxis.get_ticklines(minor=True):
        tick.set_markersize(2.5)
    for tick in ax.yaxis.get_ticklines(minor=True):
        tick.set_markersize(2.5)
