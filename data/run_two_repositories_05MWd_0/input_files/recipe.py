import numpy as np

def recipe(spent_fuel, spent_rep_fuel, burnup):
    spent_fuel_comp = np.load(spent_fuel, allow_pickle=True).item()
    spent_fuel_composition = spent_fuel_comp[burnup]
    del spent_fuel_composition['Burnup']
    del spent_fuel_composition['Irr_Time']
    
    spent_rep_fuel_comp = np.load(spent_rep_fuel, 
                                  allow_pickle=True).item()
    spent_rep_fuel_composition = spent_rep_fuel_comp[burnup]
    del spent_rep_fuel_composition['Burnup']
    del spent_rep_fuel_composition['Irr_Time']
    
    d = {"recipe": [
          {
            "name": "NaturalURecipe",
            "basis": "atom",
            "nuclide": [
              {"id": "U234", "comp": 0.0054},
              {"id": "U235", "comp": 0.7204},
              {"id": "U238", "comp": 99.2742}
            ]
          },
          {
            "name": "SRSFuelRecipe",
            "basis": "atom",
            "nuclide": [
              {"id": "U235", "comp": 0.011},
              {"id": "U238", "comp": 0.989}
            ]
          },
          {
            "name": "SpentFuelRecipe",
            "basis": "mass",
            "nuclide": [
              {"id": isotope, "comp": concentration}
                for isotope, concentration in spent_fuel_composition.items()
            ]
          },
          {
            "name": "SpentReprocessedFuelRecipe",
            "basis": "mass",
            "nuclide": [
              {"id": isotope, "comp": concentration}
                for isotope, concentration in spent_rep_fuel_composition.items()
            ]
          },
          {
            "name": "WeapongradeURecipe",
            "basis": "atom",
            "nuclide": [
              {"id": "U235", "comp": 0.9},
              {"id": "U238", "comp": 0.1}
            ]
          }
        ]}

    return d
