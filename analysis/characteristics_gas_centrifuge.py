#!/usr/bin/env python3

import matplotlib
matplotlib.use('pgf')
import matplotlib.pyplot as plt
import numpy as np

plt.style.use('seaborn')
import plotsettings as ps
plt.rcParams.update(ps.tex_fonts())


# P1 type centrifuge
CUT = 0.5
D_RHO = 2.2e-5  # kg m^-1 s^-1
DELTA_M = 3e-3  # kg mol^-1
GAS_CONSTANT = 8.314472  # m^2 kg s^-2 K^-1 mol^-1
HEIGHT = 1.8  # m
M_U238 = 238e-3  # kg mol^-1
M_UF6 = 352e-3  # kg mol^-1
MASS_RATIO = M_U238 / M_UF6
R1_R2 = 0.534  # r_1 / r_2
R2_A = 0.96  # Should lie in the range 0.96 to 0.99
S_TO_YR = 3600 * 24 * 365.25  # seconds to years
TEMPERATURE = 320  # K
VELOCITY = 320  # m s^-1

'''
# P2 type centrifuge
CUT = 0.5
D_RHO = 2.2e-5  # kg m^-1 s^-1
DELTA_M = 3e-3  # kg mol^-1
GAS_CONSTANT = 8.314472  # m^2 kg s^-2 K^-1 mol^-1
HEIGHT = 1.  # m
M_U238 = 238e-3  # kg mol^-1
M_UF6 = 352e-3  # kg mol^-1
MASS_RATIO = M_U238 / M_UF6
R1_R2 = 0.746  # r_1 / r_2
R2_A = 0.96  # Should lie in the range 0.96 to 0.99
S_TO_YR = 3600 * 24 * 365.25  # seconds to years
TEMPERATURE = 320  # K
VELOCITY = 485  # m s^-1
'''

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
def main():
    plot()
    return

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
def a_factors(feed, countercurrent_to_feed, cut):
    factor = 2 * np.pi * D_RHO * MASS_RATIO / np.log(1/R1_R2)

    a_p = (factor / feed * cut / (1+countercurrent_to_feed) 
           / (1-cut+countercurrent_to_feed))
    a_w = (factor / feed * (1-cut) / countercurrent_to_feed 
           / (1-cut+countercurrent_to_feed))
    return a_p, a_w

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
def rectifier_length(countercurrent_to_feed, cut):
    return ((1-cut) * (1+countercurrent_to_feed) * HEIGHT
            / (1 - cut +countercurrent_to_feed))

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
def separative_power(feed, countercurrent_to_feed, cut):
    a_p, a_w = a_factors(feed, countercurrent_to_feed, cut)
    z_p = rectifier_length(countercurrent_to_feed, cut)

    factor11 = 0.5 * DELTA_M * VELOCITY**2 / GAS_CONSTANT / TEMPERATURE
    factor12 = 1 - R1_R2**2
    line1 = 0.5 * feed * cut * (1-cut) * factor11**2 * factor12**2 * R2_A**4
    line2 = (1 + countercurrent_to_feed) / cut * (1-np.exp(-a_p*z_p))
    line3 = countercurrent_to_feed / (1-cut) * (1-np.exp(-a_w * (HEIGHT-z_p)))
    
    swu = line1 * (line2 + line3)**2 
    return swu

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
def separation_factor(feed, countercurrent_to_feed, cut):
    swu = separative_power(feed, countercurrent_to_feed, cut)
    alpha = 1 + (2 * swu * (1-cut) / cut / feed)**0.5
    return alpha**2

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
def plot():
    feed = np.linspace(1, 30, 50)*1e-6  # kg s^-1
    f_opt = np.array([10e-6, 4e-6]) # optimum feed rate F*, kg s^-1

    k = np.array([2.31, 3.94])
    
    fs = ps.set_size(subplots=(2, 1), fraction=0.8)
    fig, ax = plt.subplots(2, 1, figsize=fs, sharex=True)

    for i, (kk, ff) in enumerate(zip(k, f_opt)):
        countercurrent_to_feed = kk * ff / feed
        swu = separative_power(feed*MASS_RATIO, countercurrent_to_feed, CUT) * S_TO_YR
        alpha = separation_factor(feed*MASS_RATIO, countercurrent_to_feed, CUT)

        ax[0].plot(feed * 1e6, swu, label=r'$k = {}$'.format(kk),
                   linestyle=ps.linestyles(i))
        ax[1].plot(feed * 1e6, alpha, label=r'$k = {}$'.format(kk),
                   linestyle=ps.linestyles(i)) 

    ax[0].set_ylabel('Separative power [kg SWU/yr]')
    ax[0].legend()
    
    ax[1].set_xlabel('Feed rate [mg/s]')
    ax[1].set_ylabel('Stage separation factor')
    ax[1].legend()
    
    plt.tight_layout()
    plt.savefig("../plots/characteristics_swu_alpha.pdf")
    plt.close()

    return

if __name__=='__main__':
    main()
