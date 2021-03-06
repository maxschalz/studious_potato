#!/usr/bin/env python3

import matplotlib
matplotlib.use('pgf')
import matplotlib.pyplot as plt
import numpy as np

plt.style.use('seaborn-darkgrid')
import plotsettings as ps
plt.rcParams.update(ps.tex_fonts())

from tails_miso_calculator import Multi_isotope

def main():
    plot()

    return

def tails_concentration():
    """Calculate the possible enrichment given tails U234 concentrations
    """
    # Split product in two parts to have a higher granularity in the 
    # first part.
    size = 50
    xp = np.concatenate((np.linspace(1, 20, size, endpoint=False), 
                         np.linspace(20, 95, size)))
    tails = np.empty(shape=(len(xp), 2), dtype=float)

    for i, p in enumerate(xp):
        for j, u234 in enumerate((5.1e-3, 5.4e-3)):
            concentration = {'234': u234, '235': (0.720, p, 0.3)}
            m = Multi_isotope(concentration, process='centrifuge', 
                              alpha235=1.35, product=1)
            m.calculate_staging()
            tails[i,j] = m.xt[2]*100  # convert to percentage
    
    wgu_234 = tails[np.argmin(np.abs(xp-90.)), 1]
    wgu = xp[np.argmin(np.abs(xp-90.))]
    leu = xp[np.argmin(np.abs(tails[:,0]-wgu_234))]

    print(f"For x(tails,234) = {wgu_234}, the enrichment assay lies "
          + f"between {leu} and {wgu}.")
    
    return xp, tails

def plot():
    xp, tails = tails_concentration()
    tails_qty, product, u235_grade, bomb = product_qty()
 
    fig, ax = plt.subplots(figsize=ps.set_size(fraction=0.7))
    ax.plot(tails[:,1], xp, label=r'$5.4\times10^{-3}$',
               color=ps.colors(0), linestyle=ps.linestyles(0))
    ax.plot(tails[:,0], xp, label=r'$5.1\times10^{-3}$',
               color=ps.colors(1), linestyle=ps.linestyles(1))
    ax.fill_betweenx(xp, tails[50,0], tails[50,1], color='C3', 
                        alpha=0.2, label="distinction\nimpossible")
    ax.set_xlabel(r'$x_{234,T}$ [%at]')
    ax.set_ylabel(r'$x_{235,P}$ [%at]')
    ax.set_ylim(0,100)
    ax.legend(title='$x_{234,F}$ [%at]')
    
    plt.tight_layout()
    plt.savefig('../plots/tails_product_enrich.pdf')
    plt.close()


    fig, ax = plt.subplots(figsize=ps.set_size(fraction=0.8))
    for i in range(len(u235_grade)):
        ax.plot(tails_qty, product[:,i], linestyle=ps.linestyles(i), 
                   label=r'{}'.format(u235_grade[i]), color=ps.colors(i))
    
    ax.scatter(tails_qty[bomb], product[bomb, -1], marker='*', s=88,
                  color='C3', label='1 SQ U-235')
    ax.set_xlabel("tails quantity [kg]")
    ax.set_ylabel("product quantity [kg]")
    ax.set_xlim(0, tails_qty[-1])
    ax.set_yscale('log')
    ax.legend(title=r"$x_{235,P}$ [\%at]", ncol=3)
    
    plt.tight_layout()
    plt.savefig('../plots/tails_product_quantity.pdf')
    plt.close()

    return

def product_qty():
    """Determine the amount of product after finding out its assay
    """
    u235_grade = (3, 5, 20, 90)
    tails = np.linspace(100, 6500, 100)
    product = np.empty(shape=(len(tails), len(u235_grade)), dtype=float) 

    for i, t in enumerate(tails):
        for j, u235 in enumerate(u235_grade):
            concentration = {'234': 5.4e-3, '235': (0.711, u235, 0.3)}
            m = Multi_isotope(concentration, process='centrifuge', 
                              alpha235=1.3, tails=t)
            m.calculate_staging()
            product[i,j] = m.p
   
    bomb = np.argmin(np.abs(product[:,-1]-27.8))
    
    return tails, product, u235_grade, bomb
    
def product_qty_unknown_feed():
    """Determine the amount of product after finding out its assay
    """

    u235_grade = (3, 20, 90)
    tails = np.linspace(100, 6500, 100)
    product = np.empty(shape=(len(tails), len(u235_grade), 2), dtype=float) 

    for i, t in enumerate(tails):
        for j, u235 in enumerate(u235_grade):
            for k, u234 in enumerate((5.4e-3, 5.1e-3)):
                concentration = {'234': u234, '235': (0.711, u235, 0.3)}
                m = Multi_isotope(concentration, process='centrifuge', 
                                  alpha235=1.3, tails=t)
                m.calculate_staging()
                product[i,j,k] = m.p
       
    fig, ax = plt.subplots(figsize=ps.set_size())
    for i in range(len(u235_grade)):
        ax.plot(tails, product[:,i,0], label=r'{}'.format(u235_grade[i]),
                linestyle=ps.linestyles(i), color=ps.colors(i))
        ax.fill_between(tails, product[:,i,0], product[:,i,1],
                        linestyle=ps.linestyles(i), color=ps.colors(i))
    
    bomb = np.argmin(np.abs(product[:,-1,0]-27.8))
    ax.scatter(tails[bomb], product[bomb,-1,0], marker='*', s=88,
               color='C2', label='1 SQ U-235')
    ax.set_xlabel("tails quantity [kg]")
    ax.set_ylabel("product quantity [kg]")
    ax.set_xlim(0, tails[-1])
    ax.set_yscale('log')
    ax.legend(title=r"$x_{235,P}$ [\%at]")
    
    plt.tight_layout()
    plt.savefig('../plots/productqty_vs_tailsqty_feed_uncertainty.pdf')
    plt.close() 

    return
   
if __name__=='__main__':
    main()
