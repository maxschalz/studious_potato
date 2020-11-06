#!/usr/bin/env python3

import numpy as np
import re
from scipy.optimize import minimize

from analyser import Analyser
from multi_isotope_calculator import Multi_isotope

# Global variables
NUC_ID = (922320000, 922330000, 922340000, 922350000, 922360000, 922380000)
SEPARATION_FACTOR = 1.35

def main():
    burnup = ("0.5MWd", "2MWd")
    origin = ("Natural", "Reprocessed")
    
    for bu in burnup:
        for orig in origin:
            tails_qty(orig, bu)
    return

def simulation_tails(fname, uranium_type="Natural"):
    """Get the reprocessed and depleted uranium tails as dict"""
    
    if uranium_type=="Natural":
        sink = "DepletedNaturalUSink"
    elif uranium_type=="Reprocessed":
        sink = "DepletedReprocessedUSink"
    else:
        msg = "'uranium_type' has to be either 'Natural' or 'Reprocessed'."
        raise ValueError(msg)
    
    a = Analyser(fname)
    sim_end = a.query(selection='EndTime', table='Finish')[0][0]
    results = a.query(selection='NucId, Quantity',
                      table='ExplicitInventory',
                      condition='Time==? AND AgentId==?',
                      vals=(sim_end, a.names[sink]))
    comp = dict([(key, 0) for key in range(232, 239) if key!=237] )
    quantity = 0
    for key, value in results:
        key = key/10000 - 92*1000
        comp[key] = value
        quantity += value

    for key, value in comp.items():
        comp[key] = value / quantity
    
    return comp, quantity

def enrichment_feed_and_tails(origin, burnup):
    """Prepare data: feed used and expected tails composition"""
    if origin=="Natural":
        #Mass fractions of natU enrichment tails taken from the Cyclus 
        # output file
        natU_comp = {'234': 0.0054, '235': (0.7204, 2, 0.3)}
        seu_tails = np.array([0., 0., 1.5440247618063e-05, 
                              0.00290322658192604, 0., 
                              0.997081333170456])
        heu_tails = np.array([0., 0., 1.27218682709261e-05, 
                              0.00285479562964945, 0., 
                              0.99713248250208])
        return natU_comp, (seu_tails, heu_tails)
    elif origin=="Reprocessed":
        # Load, filter and format feed data
        data = np.load("../data/SERPENT_outputs_NatU_percentages.npy").item()

        feed_composition = {}
        normalisation = 0
        for iso in [i for i in range(234, 239) if i!=237]:
            value = data[burnup][f"U{iso}"]
            feed_composition[str(iso)] = value * 100
            normalisation += value

        for key, val in feed_composition.items():
            feed_composition[key] = val/normalisation
        
        feed_composition['232'] = 0.
        feed_composition['233'] = 0
        # The U238 content is calculated by the enrichment module
        del feed_composition['238']
       
        # Get SEU and HEU tails
        if burnup=="0.5MWd":
            seu_tails = np.array([0., 0., 1.35406410557832e-05, 
                                  0.00269133511129306, 4.13592084547905e-05,
                                  0.997253765039196])
        elif burnup=="2MWd":
            seu_tails = np.array([0., 0., 1.56662456546925e-05, 
                                  0.00269248329581373, 0.000163308471630726,
                                  0.997128541986901])
            
        else:
            raise ValueError("'burnup' has to be '0.5MWd' or '2MWd'")
        concentration = feed_composition
        concentration['235'] = (feed_composition['235'], 90., 0.3)
        m = Multi_isotope(concentration, feed=1, process='centrifuge', 
                          alpha235=SEPARATION_FACTOR, downblend=True)
        m.calculate_staging()
        heu_tails = m.xt
        
        return feed_composition, (seu_tails, heu_tails)
    else:
        raise ValueError("'origin' has to be 'Natural' or 'Reprocessed'")

def mix_compositions(comp1, comp2, mass1, mass2):
    return (mass1*comp1 + mass2*comp2) / (mass1+mass2)

def mixing_ratios(sim_tails_comp, tails_comp):
    """Calculate how much of comp1 is added to comp2 using mix_comp
    
    Here, the mixing with the following compositions is calculated:
    mix_comp = (a*comp1 + b*comp2) / (a+b)
    b is set to 1 such that this function calculates how much of
    comp1 is added to comp2 per unit of comp2. In other words, a is
    given in units of comp2.
    """
    # Assure correct formatting
    sim_tails_comp = np.array(list(sim_tails_comp.values()))
        
    # special case: mix comp contains no comp2: return a large number
    if np.all(sim_tails_comp - tails_comp[0] < 1e-10):
        print(f"Only SEU tails, no HEU tails produced!")
        return 1e200
    
    mass_ratio = ((tails_comp[1]-sim_tails_comp) 
                  / (sim_tails_comp-tails_comp[0]))

    if np.std(mass_ratio[~np.isnan(mass_ratio)]) > 1e-10:
        print()
        msg = (f"Values differ from each other!\n"
               + f"mass_ratio:\n{mass_ratio}\n"
               + f"Composition 1:\n{tails_comp[0]}\n"
               + f"Composition 2:\n{tails_comp[1]}\n"
               + f"Mixed final composition:\n{sim_tails_comp}")
        raise RuntimeError(msg)
    
    # Remove possible nans from isotopes
    mass_ratio = np.mean(mass_ratio[~np.isnan(mass_ratio)])
        
    return mass_ratio

def tails_per_product_qty(concentrations, enrichment_level):
    m = Multi_isotope(concentrations, max_swu=np.inf, feed=np.inf, 
                      product=1, downblend=True, process='centrifuge',
                      alpha235=SEPARATION_FACTOR)
    m.set_product_enrichment(enrichment_level)
    m.calculate_staging()
    tails = m.t
    product = m.p
    if abs(product-1) > 1e-10:
        raise RuntimeError("Something fishy going on here")
    
    return tails / product

def tails_qty(origin, burnup):
    """Calculate the amount of HEU and SEU produced"""
    
    print(f"\n{origin} uranium, burnup of {burnup}")
    # Get tails composition in depleted U sink from simulation
    fname_burnup = re.sub("\.", "", burnup)
    fname = (f"../data/run_two_repositories_{fname_burnup}_0/"
             + f"run_two_repositories_{fname_burnup}.sqlite")
    sim_tails_comp, sim_tails_qty = simulation_tails(fname, 
                                                     uranium_type=origin)

    # Get feed and predicted tails compositions
    feed_comp, tails_comp = enrichment_feed_and_tails(origin, burnup)
    
    seu_per_heu_tails = mixing_ratios(sim_tails_comp, tails_comp)    
    seu_tails_qty = (sim_tails_qty * seu_per_heu_tails 
                     / (1.+seu_per_heu_tails))
    heu_tails_qty = sim_tails_qty / (1.+seu_per_heu_tails)
    
    print(f"Total qty:    {sim_tails_qty:9.0f} kg\n"
          + f"SEU tails:    {seu_tails_qty:9.0f} kg\n"
          + f"HEU tails:    {heu_tails_qty:9.0f} kg\n")
    
    enrichment_lvl = (1.1, 90.)
    label = ("SEU", "HEU")
    tails = (seu_tails_qty, heu_tails_qty)
    
    for xp, name, tail in zip(enrichment_lvl, label, tails):
        t_per_p = tails_per_product_qty(feed_comp, xp)
        product = tail / t_per_p
        print(f"Produced {name}: {product:9.1f} kg")
    
    return

if __name__=="__main__":
    main()
