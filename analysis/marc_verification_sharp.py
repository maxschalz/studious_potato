#!/usr/bin/env python3

import matplotlib
matplotlib.use('pgf')
import matplotlib.pyplot as plt
import numpy as np

from multi_isotope_calculator import Multi_isotope

import plotsettings as ps
plt.style.use('seaborn-darkgrid')
plt.rcParams.update(ps.tex_fonts())


def main():
    plot()
    #figure5()

def figure1():
    """Compare data to Sharp paper (tails U234 vs product U235)"""
    data = np.genfromtxt("../data/sharp_fig1.csv", delimiter=",")
    data = data[np.argsort(data[:,0])]

    composition = {'234': 5.5e-3, '235': (0.72, 3, 0.2)}
    calculator = Multi_isotope(composition, feed=1, process='diffusion',
                               downblend=False)
    results = np.empty(shape=data.shape, dtype=float)
    for i, xp in enumerate(data[:,0]):
        calculator.set_product_enrichment(xp*100)
        calculator.calculate_staging()
        results[i,0] = calculator.xp[3]
        results[i,1] = calculator.xt[2]

    data *= 100
    results *= 100
    pulls = 100 * (data[:,1]-results[:,1]) / data[:,1]    
    
    ylims = (1e299, 0)
    for values in (data, results):
        ylims = (min(ylims[0], min(values[:,1])),
                 max(ylims[1], max(values[:,1]))) 
    
    return data, results, pulls

def figure5():
    """Compare data to Sharp paper (tails qty vs product qty)"""
    sharp = np.genfromtxt("../data/sharp_fig5.csv", delimiter=",")
    sharp = sharp[np.argsort(sharp[:,0])]
    
    calc = Multi_isotope({'235': (0.711, 5, 0.2)}, max_swu=15000,
                         process='diffusion', downblend=False)

    results = np.empty(shape=sharp.shape, dtype=float)
    for i, xp in enumerate(sharp[:,0]):
        calc.set_product_enrichment(xp*100)
        calc.calculate_staging()
        
        results[i,0] = calc.xp[3] * 100
        results[i,1] = calc.t
   
    sharp[:,0] *= 100
    pulls = 100 * (sharp[:,1]-results[:,1]) / sharp[:,1]

    return sharp, results, pulls 

def plot():
    fig1 = figure1()
    fig5 = figure5()

    figsize = ps.set_size(subplots=(2,2))
    fig, ax = plt.subplots(figsize=figsize, nrows=2, ncols=2)
    
    plt.rcParams.update({'lines.markersize': 4})

    for i, (data, result, pulls) in enumerate((fig1, fig5)):
        ax[0,i].plot(result[:,0], result[:,1], color=ps.colors(0),
                     label="MARC algorithm", zorder=2, linewidth=1)
        ax[0,i].scatter(data[::3,0], data[::3,1], marker="x",
                        color=ps.colors(1), label="Sharp 2013", zorder=3)
        ax[1,i].scatter(data[:,0], pulls, s=1, zorder=2)
        
        ax[0,i].legend()
        ax[0,i].set_xlim(0, 100)
        ax[1,i].set_xlim(0, 100)
        ax[1,i].set_xlabel(r"$x_{235,P}$ [\%at]")
        ax[1,i].axhline(0, color="C3", zorder=1, linewidth=1)
    
    ax[0,1].ticklabel_format(axis="y", style="sci", scilimits=(-2,2))
    ax[0,0].set_ylabel(r"$x_{234,T}$ [\%at]") 
    ax[1,0].set_ylabel(r"relative difference [%]")
    ax[0,1].set_ylabel(r"$T$ [kg/yr]") 
    ax[1,1].set_ylabel(r"relative difference [%]")

    plt.tight_layout()
    plt.savefig("../plots/checks_marc_sharp1.pdf")
    plt.close()
    

    return


if __name__=='__main__':
    main()
