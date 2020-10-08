#!/usr/bin/env python3

import numpy as np
import os
import sys
sys.path.append('../../enrichment/src')

from multi_isotope_calculator import Multi_isotope

DATA_PATH = '../data/'
CORE_MASS = 110820  # in kg
REFUELLING_TIME = 6  # in days
SEPARATION_EFFICIENCY = 0.97 
SIM_DUR = 729  # in days

def main():
    spent_natU_fname = os.path.join(DATA_PATH,
            "SERPENT_outputs_NatU_percentages.npy")

    get_expectations(spent_natU_fname, 22, '0.5MWd', True)
    get_expectations(spent_natU_fname, 88, '2MWd', True)

    return

def get_expectations(fname, irradiation_time, burnup, verbose=True):
    """Run all the calculations for one simulation""" 
    
    print(f"Expected values for a cycle with burnup {burnup} and "
          + f"irradiation time {irradiation_time}:\n")

    reactor_cycles(irradiation_time, verbose=True)
    natU_to_repU_cycles(fname, irradiation_time, burnup, verbose)
    expected_plutonium(burnup, irradiation_time)
    expected_heu(fname, irradiation_time, burnup)
    spent_reprocessed_uranium(fname, burnup, irradiation_time, verbose)
    print("\n\n")
    
    return

def expected_heu(fname, irradiation_time, burnup):
    """Get the expected amount of weapongrade U produced"""

    # Get fraction of times where natU is used as reactor fuel
    n_natU, n_repU = natU_to_repU_cycles(fname, irradiation_time, burnup,
                                         verbose=False)

    # + 1 because of the stored fuel assembly and + 1 for the incomplete
    # reactor cycle at the end of the simulation.
    time_reactor_enrichment = (n_natU + 2) * enrichment_reactorgrade()
    time_heu_enrichment = SIM_DUR - time_reactor_enrichment

    m = Multi_isotope({'234': 0.0054, '235': (0.7204, 90., 0.3)},
                      feed=10000, alpha235=1.35, process='centrifuge',
                      downblend=True)
    m.calculate_staging()
    heu_per_cycle = m.p

    total_heu = heu_per_cycle * time_heu_enrichment

    print((f'Total weapongrade U: {total_heu:.1f} kg, using an irradiation'
          + f' time of {irradiation_time}'))
    
    return total_heu

def enrichment_reactorgrade(verbose=False):
    """Get the (non-int) timesteps needed to produce one SRS core
    
    This function returns the number of enrichment cycles needed to enrich 
    NatU to 1.1% making it usable in the Savannah River Site reaction. One
    full reactor core contains 110820 kg SEU and we assume that in one step
    10'000 kg of uranium are used as feed.

    While timesteps typically are integers, this is not the case here. 
    Using floats has the advantage that the enrichment in the last step is
    reflected better, as the facility retains some capacity that can 
    subsequently be used to enrich natural uranium to HEU.
    """
    m = Multi_isotope({'234': 0.0054, '235': (0.7204, 1.1, 0.3)},
                      feed=10000, alpha235=1.35, process='centrifuge',
                      downblend=True)
    m.calculate_staging()
    product = m.p
    
    n_steps = CORE_MASS / product
    if verbose:
        print(f'{n_steps} enrichment cycles needed.')

    return n_steps

def spent_reprocessed_uranium(fname, burnup, irradiation_time, 
                              verbose=False):
    if burnup=='0.5MWd':
        fname = os.path.join(DATA_PATH, 'SERPENT_outputs_RepU_05MWd_percentages.npy')

    elif burnup=='2MWd':
        fname = os.path.join(DATA_PATH, 'SERPENT_outputs_RepU_2MWd_percentages.npy')
    else:
        raise ValueError("'burnup' has to be either '0.5MWd' or '2MWd'")
    
    data = np.load(fname, allow_pickle=True).item()
    data = data[burnup]
    
    uranium_content = 0
    for key, val in data.items():
        if key in [f'U{iso}' for iso in range(232, 239)]:
            uranium_content += val

    spent_batch = CORE_MASS * SEPARATION_EFFICIENCY * uranium_content

    n_reactor_cycles = natU_to_repU_cycles(fname, irradiation_time, 
                                           burnup)[1]
    spent_reprocessed = n_reactor_cycles * spent_batch

    if verbose:
        print(f'{spent_reprocessed/1000:.1f} t of spent uranium in storage')

    return spent_reprocessed

def reactor_cycles(irradiation_time, verbose=False):
    """Get the (int) number of expected reactor cycles in one simulation"""

    # Note the integer division!
    n_cycles = ((SIM_DUR - enrichment_reactorgrade()) 
                // (irradiation_time + REFUELLING_TIME))
    
    if verbose:
        print(f'Irradiation time of {irradiation_time} yields {n_cycles} '
               + 'cycles')

    return n_cycles

def natU_to_repU_cycles(fname, irradiation_time, burnup, verbose=False):
    """Get the number of reactor cyclues using repU and natU"""
    data = np.load(fname, allow_pickle=True).item()
    
    data = data[burnup]

    spentU_composition = {}
    uranium_content = 0
    # Get uranium content and isotopic composition of uranium
    for key, val in data.items():
        if key in [f'U{iso}' for iso in range(232, 239) if iso!=237]:
            spentU_composition[key[1:]] = val  # remove the 'U' 
            uranium_content += val
    
    # Normalise spent uranium's isotopic composition to 100 (percent)
    for key, val in spentU_composition.items():
        spentU_composition[key] = 100 * val / uranium_content

    total_cycles = reactor_cycles(irradiation_time, False)
    spent_batch = CORE_MASS * SEPARATION_EFFICIENCY * uranium_content
    
    spentU_composition['235'] = (spentU_composition['235'], 1.1, 0.3)
    del spentU_composition['238']
    m = Multi_isotope(spentU_composition, feed=spent_batch, alpha235=1.35,
                      process='centrifuge', downblend=True)
    m.calculate_staging()
    repU_batch_mass = m.p
    
    cycle = 0  # timestep in the form of reactorcycles
    repU_storage = 0
    n_natU = 0
    n_repU = 0
    while cycle < total_cycles:
        if repU_storage < CORE_MASS:
            n_natU += 1
            repU_storage += repU_batch_mass
        else:
            n_repU += 1
            repU_storage -= CORE_MASS
        cycle += 1

    if verbose:
        print(f"Out of {total_cycles} cycles, {n_natU} used fresh fuel and"
              + f" {n_repU} used reprocessed fuel.")
        if repU_storage < CORE_MASS:
            print(f"Last, incomplete cycle uses natU")
        else:
            print(f"Last, incomplete cycle uses natU")

    return (n_natU, n_repU)

def expected_plutonium(burnup, irradiation_time):
    """This is ugly coding don't look at it"""

    data = [] 
    data.append(get_plutonium(
        os.path.join(DATA_PATH, 'SERPENT_outputs_NatU_percentages.npy'),
        burnup))
    if burnup=='0.5MWd':
        pu = get_plutonium(os.path.join(DATA_PATH, 
                'SERPENT_outputs_RepU_05MWd_percentages.npy'), burnup)
        data.append(pu)
    elif burnup=='2MWd':
        pu = get_plutonium(os.path.join(DATA_PATH, 
                'SERPENT_outputs_RepU_2MWd_percentages.npy'), burnup)
        data.append(pu)
    else:
        raise ValueError("'burnup' has to be either '0.5MWd' or '2MWd'")
    
    plutonium = np.array(data) 
    plutonium *= (reactor_cycles(irradiation_time) * CORE_MASS 
                  * SEPARATION_EFFICIENCY)
    mean = np.mean(plutonium)
    std = np.std(plutonium, ddof=1)
    print(f"Pu for {burnup}: {mean:.1f} +- {std:.1f}")

    return

def get_plutonium(fname, burnup):
    """Load the composition of spent fuel and filter it (e.g., only U)"""
    data = np.load(fname, allow_pickle=True).item()
    data = data[burnup]
    
    pu = 0
    for isotope, value in data.items():
        if isotope in ('Pu239', 'Pu240', 'Pu241'):
            pu += value

    pu += (1 - 2.**(-1./2.356)) * data['Np239']

    return pu

if __name__=="__main__":
    main()
