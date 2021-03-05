import os
import sys
sys.path.append("./")

import archetypes
import commodity
import control
import facility
import institution
import recipe
import region

def simulation():
    data_path = "../.."

    burnup = '0.5MWd'
    #burnup = '2MWd'

    spent_fuel_fname = "SERPENT_outputs_NatU_percentages.npy"
    spent_fuel_fname = os.path.join(data_path, spent_fuel_fname)
    
    if burnup == '0.5MWd':
        spent_rep_fuel_fname = "SERPENT_outputs_RepU_05MWd_percentages.npy"
    elif burnup == '2MWd':
        spent_rep_fuel_fname = "SERPENT_outputs_RepU_2MWd_percentages.npy"
    spent_rep_fuel_fname = os.path.join(data_path, spent_rep_fuel_fname)
    
    arch = archetypes.archetypes()
    commod = commodity.commodity()
    ctrl = control.control()
    fac = facility.facility(burnup)
    inst = institution.institution(burnup)  # not needed but it should be
                                            # called during the consistency
                                            # check
    recipes = recipe.recipe(spent_fuel=spent_fuel_fname,
                            spent_rep_fuel=spent_rep_fuel_fname, 
                            burnup=burnup)
    reg = region.region(burnup)
    
    return {"simulation": {**arch, **commod, **ctrl, **fac, **recipes, 
                           **reg}}

if __name__=='__main__':
    simulation()
