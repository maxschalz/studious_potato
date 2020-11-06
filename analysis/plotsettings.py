import matplotlib.pyplot as plt

def tex_fonts():
    return {
        # Use LaTeX to write all text
        "pgf.texsystem": "pdflatex",
        "pgf.preamble":  [
            r"\usepackage[utf8x]{inputenc}",    # use utf8 fonts 
            r"\usepackage[T1]{fontenc}",        # plots will be generated
            r"\usepackage[detect-all,locale=DE]{siunitx}",
            r"\usepackage[version=4]{mhchem}",
            r"\usepackage{isotope}",
        ],
        "text.usetex": True,
        "font.family": "serif",
        # Use 10pt font in plots, to match 10pt font in document
        "axes.labelsize": 10,
        "font.size": 10,
        # Make the legend/label fonts a little smaller
        "legend.fontsize": 8,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,

        "legend.fancybox": True,
        "legend.frameon": True,
        "legend.framealpha": 1,
        
        "axes.formatter.limits": (-3, 5),
        "axes.formatter.use_mathtext": True
    }

def set_size(width='book', fraction=1, subplots=(1, 1), higher=False):
    """Set figure dimensions to avoid scaling in LaTeX.

    Parameters
    ----------
    width: float or string referring to a documentclass
    Document textwidth or columnwidth in pts
    fraction: float, optional
    Fraction of the width which you wish the figure to occupy

    Returns
    -------
    fig_dim: tuple
    Dimensions of figure in inches
    """

    if width == 'book':
       width_pt = 345.0
    else:
        width_pt = width

    # Width of figure (in pts)
    fig_width_pt = width_pt * fraction

    # Convert from pt to inches
    inches_per_pt = 1 / 72.27

    # Golden ratio to set aesthetic figure height
    # https://disq.us/p/2940ij3
    golden_ratio = (5**.5 - 1) / 2

    # Figure width in inches
    fig_width_in = fig_width_pt * inches_per_pt
    # Figure height in inches

    if higher:
        fig_height_in = (fig_width_in * golden_ratio 
                         * (subplots[0] * 1.5 / subplots[1]))
    else:
        fig_height_in = (fig_width_in * golden_ratio 
                         * (subplots[0] / subplots[1]))
    fig_dim = (fig_width_in, fig_height_in)

    return fig_dim

def colors(i):
    """Return a pyplot default color
    """
    return f'C{i}'

def linestyles(i):
    """Return a linestyle
    """
    ls = ('dashed', 'dashdot', 'dotted', (0, (3, 1, 1, 1, 1, 1)),
          (0, (3, 5, 1, 5)))

    return ls[i % len(ls)]

def markers(i, markersize=4):
    """Return a marker style
    """
    plt.rcParams.update({'lines.markersize': markersize})
    mk = ('o', 'P', '^', 'X', 'p')

    return mk[i % len(mk)]
