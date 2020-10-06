#!/usr/bin/env python3

import numpy as np
import os
import sys
sys.path.append('../../enrichment/src')

from multi_isotope_calculator import Multi_isotope

CORE_MASS = 110820  # in kg
REFUELLING_TIME = 6  # in days
SEPARATION_EFFICIENCY = 0.97 
SIM_DUR = 730  # in days

def main():
    #expected_plutonium()
    #expected_heu(22)
    #expected_heu(88)
    spent_reprocessed_uranium(22, verbose=False)
    spent_reprocessed_uranium(88, verbose=False)

    return

def expected_heu(irradiation_time):
    """Get the expected amount of weapongrade U produced"""
    
    # + 1 because of the stored fuel assembly and + 1 for the incomplete
    # reactor cycle at the end of the simulation.
    n_reactor_cycles = reactor_cycles(irradiation_time) + 2
    fraction_natU_to_repU = 5 / 7
    time_reactor_enrichment = (n_reactor_cycles * enrichment_reactorgrade()
                               * fraction_natU_to_repU)
    time_heu_enrichment = SIM_DUR - time_reactor_enrichment
    
    m = Multi_isotope({'234': 0.0054, '235': (0.7204, 90., 0.3)},
                      feed=10000, alpha235=1.35, process='centrifuge',
                      downblend=True)
    m.calculate_staging()
    heu_per_cycle = m.p
    total_heu = heu_per_cycle * time_heu_enrichment

    print(time_heu_enrichment)
    print(f'Weapongrade U per enrichment step: {heu_per_cycle:.3f} kg')
    print((f'Total weapongrade U: {total_heu:.1f} kg, using an irradiation'
          + f' time of {irradiation_time}'))
    
    return total_heu

def enrichment_reactorgrade(verbose=False):
    """Get the number of enrichment cycles needed to produce one SRS core
    
    This function returns the number of enrichment cycles needed to enrich 
    NatU to 1.1% making it usable in the Savannah River Site reaction. One
    full reactor core contains 110820 kg SEU and we assume that in one step
    10'000 kg of uranium are used as feed.
    """
    m = Multi_isotope({'234': 0.0054, '235': (0.7204, 1.1, 0.3)},
                      feed=10000, alpha235=1.35, process='centrifuge',
                      downblend=True)
    m.calculate_staging()
    product = m.p

    n_steps = int(CORE_MASS / product) + 1

    if verbose:
        print(f'{n_steps} enrichment cycles needed.')

    return n_steps

def spent_reprocessed_uranium(irradiation_time, verbose=False):
    # Uranium content taken from SERPENT output file
    uranium_content = (0.988404066305649 + 0.0102986527882519 
                       + 0.000198774018864621 + 8.88584554656443e-05
                       + 1.07743552265805e-06 + 5.38515760293621e-07)
    spent_batch = CORE_MASS * SEPARATION_EFFICIENCY * uranium_content
    if verbose:
        print(f'{spent_batch:.1f} kg of stored uranium per core')

    n_reactor_cycles = reactor_cycles(irradiation_time)
    # Minus 1 to account for one missing cycle
    spent_reprocessed = int(0.5 * n_reactor_cycles - 1) * spent_batch
    print(f'{spent_reprocessed:.1f} kg of spent uranium in storage')

    return spent_reprocessed

def reactor_cycles(irradiation_time, verbose=False):
    """Get the (int) number of expected reactor cycles in one simulation"""

    # Note the integer division!
    n_cycles = ((SIM_DUR - enrichment_reactorgrade()) 
                // (irradiation_time + REFUELLING_TIME))
    
    if verbose:
        print((f'Irradiation time of {irradiation_time} yiels {n_cycles} '
               + 'cycles'))

    return n_cycles

def expected_plutonium():
    """This is ugly coding don't look at it"""

    path = '../data/'
    
    pu = [] 
    pu.extend(get_composition(
        os.path.join(path, 'SERPENT_outputs_NatU_percentages.npy')))
    pu.extend(get_composition(
        os.path.join(path, 'SERPENT_outputs_RepU_05MWd_percentages.npy'),
        bu='0.5MWd'))
    pu.extend(get_composition(
        os.path.join(path, 'SERPENT_outputs_RepU_2MWd_percentages.npy'),
        bu='2MWd'))
    
    n_enrichment_steps = enrichment_reactorgrade()
    get_plutonium = lambda pu, n_cycles: (pu * n_cycles * CORE_MASS
                                          * SEPARATION_EFFICIENCY)
 
    low_bu_natU = get_plutonium(pu[0], reactor_cycles(22))
    low_bu_repU = get_plutonium(pu[2], reactor_cycles(22))
    
    high_bu_natU = get_plutonium(pu[1], reactor_cycles(88))
    high_bu_repU = get_plutonium(pu[3], reactor_cycles(88))
                                        
    print('Pu for 0.5MWd: {} -- {}'.format(low_bu_repU, low_bu_natU))
    print('Pu for 2.0MWd: {} -- {}'.format(high_bu_repU, high_bu_natU))

    return

def get_composition(fname, bu=None):
    """Load the composition of spent fuel and filter it (e.g., only U)"""
    
    fuel = fname.split('_')[2]
    burnup_fname = fname.split('_')[-2]

    print(f"Extracting data for {fuel}, {burnup_fname}")
    data = np.load(fname, allow_pickle=True).item()
   
    bu = data.keys() if bu is None else [bu]
    pu = []
    for burnup in bu:
        comp = {'waste': 0, 'U234': 0, 'U235': 0, 'U236': 0, 'U238': 0,
                'Pu239': 0, 'Pu240': 0, 'Pu241': 0, 'Np239': 0}
        print(f"Extracting data for burnup of {burnup}")
        del data[burnup]['Burnup']
        del data[burnup]['Irr_Time']
        
        for isotope, value in data[burnup].items():
            if isotope not in comp.keys():
                comp['waste'] += value
            else:
                comp[isotope] += value

        #print(comp)
        pu_comp = sum([data[burnup][key] 
                          for key in ('Pu239', 'Pu240', 'Pu241')])
        pu_comp += (1 - 2.**(-1./2.356)) * data[burnup]['Np239']
        #pu_comp += 2.**(-1./2.356) * data[burnup]['U237']
        print('Plutonium: {} %\n'.format(pu_comp*100))

        pu.append(pu_comp)
    
    return pu

if __name__=="__main__":
    main()
