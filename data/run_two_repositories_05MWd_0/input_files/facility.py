import sys
sys.path.append(
    "/Users/test/Uni/Masterarbeit/thesis_simulations/input/test_transfers/"
)

def facility(burnup):
    if burnup == '0.5MWd':
        cycle_time = 22  # days
    elif burnup == '2MWd':
        cycle_time = 88  # days
    else:
        raise ValueError("Burnup value has to be '0.5MWd' or '2MWd'!")

    refuel_time = 6  # days

    # The masses below are all indicated in kg and the concept is 
    # constructed such that a full core always weighs `reactor_mass` kg and
    # is completely changed at the end of a cycle.
    reactor_mass = 110.82e3  # mass of the full core
    
    n_assem_batch = 1   # number of assemblies per batch, number of 
                        # assemblies discharged from the core fully burned 
                        # each cycle.
    assembly_size = reactor_mass / n_assem_batch  # mass of a single 
                                                  # assembly
    n_assem_core = n_assem_batch  # number of assemblies in a full core
    n_assem_fresh = 0#n_assem_batch  # keep one fresh batch at hand to ensure 
                                   # continuous operation of the reactor
    n_assem_spent = n_assem_batch  # The NFC is constructed such that spent
                                   # assemblies always get transported away
                                   # so this should never be a problem.
    power_thermal = 2496  # in MW
    power_electric = power_thermal / 3.  # in MW

    gamma_235 = 1.35
    swu_per_factory = 1e299
    source_throughput = 10000
    
    separation_efficiency = 0.97  # Taken from 'The nuclear fuel cycle',
                                  # ed. by P. D. Wilson, Oxford University
                                  # Press (1996): p. 138

    d = {"facility": [
          {
            "name": "NaturalUSource",
            "config": {"Source": {
              "outcommod": "NaturalU",
              "outrecipe": "NaturalURecipe",
              "throughput": source_throughput
            }}
          },
          {
            "name": "NaturalUEnricher",
            "config": {"MIsoEnrich": {
              "feed_commod": "NaturalU",
              "feed_recipe": "NaturalURecipe",
              "product_commod": "EnrichedU",
              "tails_commod": "DepletedNaturalU",
              "tails_assay": 0.003,
              "order_prefs": False,
              "initial_feed": 0,
              "max_feed_inventory": 1e299,  #reactor_mass*2.4,
              "gamma_235": gamma_235,
              "swu_capacity": swu_per_factory,
              "swu_capacity_vals": {"val": [swu_per_factory]},
              "swu_capacity_times": {"val": [0]},
              "use_downblending": True
            }}
          },
          {
            "name": "NatUFuelStorage",
            "config": {"Storage": {
              "in_commods": {"val": ["EnrichedU"]},
              "out_commods": {"val": ["SRSFuel"]},
              "in_recipe": "SRSFuelRecipe",
              "residence_time": 0,
              "throughput": 1e299,
              "max_inv_size": reactor_mass
            }}
          },
          {
            "name": "ReprocessedEnricher",
            "config": {"MIsoEnrich": {
              "feed_commod": "ReprocessedU",
              "feed_recipe": "NaturalURecipe",
              "product_commod": "ReprocessedFuel",
              "tails_commod": "DepletedReprocessedU",
              "tails_assay": 0.003,
              "max_enrich": 0.1,
              "order_prefs": False,  
              "initial_feed": 0,
              "max_feed_inventory": 1e299,
              "gamma_235": gamma_235,
              "swu_capacity": swu_per_factory,
              "swu_capacity_vals": {"val": [swu_per_factory]},
              "swu_capacity_times": {"val": [0]},
              "use_downblending": True
            }}
          },
          {
            "name": "RepUFuelStorage",
            "config": {"Storage": {
              "in_commods": {"val": ["ReprocessedFuel"]},
              "out_commods": {"val": ["SRSReprocessedFuel"]},
              "in_recipe": "SRSFuelRecipe",
              "residence_time": 0,
              "throughput": 1e299,
              "max_inv_size": 1e299,
            }}
          },
          {
            "name": "SRS",
            "config": {"Reactor": {
              "fuel_incommods": {"val": ["SRSReprocessedFuel", "SRSFuel"]},#["ReprocessedFuel", "EnrichedU"]},
              "fuel_inrecipes": {"val": ["SRSFuelRecipe", "SRSFuelRecipe"]},
              "fuel_prefs": {"val": [5, 1]},
              "fuel_outcommods": {"val": ["SpentReprocessedFuel", "SpentFuel"]},
              "fuel_outrecipes": {"val": ["SpentReprocessedFuelRecipe",
                                          "SpentFuelRecipe"]},
              "assem_size": assembly_size,
              "n_assem_batch": n_assem_batch,
              "n_assem_core": n_assem_core,
              "n_assem_fresh": n_assem_fresh,
              "n_assem_spent": n_assem_spent,
              "cycle_time": cycle_time,
              "refuel_time": refuel_time,
              "power_cap": power_electric  # MW_e
            }}
          },
          {
            "name": "Recycler",
            "config": {"Separations": {
              "feed_commods": {"val": ["SpentFuel"]},
              "feedbuf_size": 1e8,
              "leftover_commod": "NuclearWaste",
              "streams": {
                "item": [
                  {
                    "commod": "WeapongradePu",
                    "info": {
                      "buf_size": 1e299,
                      "efficiencies": {
                        "item": [{
                          "comp": "94000",
                          "eff": separation_efficiency
                        }]
                      }
                    }
                  },
                  {
                    "commod": "ReprocessedU",
                    "info": {
                      "buf_size": 1e299,
                      "efficiencies": {
                        "item": [{
                          "comp": "92232",
                          "eff": separation_efficiency
                        },
                        {
                          "comp": "92233",
                          "eff": separation_efficiency
                        },
                        {
                          "comp": "92234",
                          "eff": separation_efficiency
                        },
                        {
                          "comp": "92235",
                          "eff": separation_efficiency
                        },
                        {
                          "comp": "92236",
                          "eff": separation_efficiency
                        },
                        {
                          "comp": "92238",
                          "eff": separation_efficiency
                        }]
                      }
                    }
                  }
                ]
              }
            }}
          },
          {
            "name": "Plutoniumseparator",
            "config": {"Separations": {
              "feed_commods": {"val": ["SpentReprocessedFuel"]},
              "feedbuf_size": 1e8,
              "leftover_commod": "NuclearWaste",
              "streams": {
                "item": [
                  {
                    "commod": "WeapongradePu",
                    "info": {
                      "buf_size": 1e299,
                      "efficiencies": {
                        "item": [{
                          "comp": "94000",
                          "eff": separation_efficiency
                        }]
                      }
                    }
                  },
                  {
                    "commod": "SpentU",
                    "info": {
                      "buf_size": 1e299,
                      "efficiencies": {
                        "item": [{
                          "comp": "92232",
                          "eff": separation_efficiency
                        },
                        {
                          "comp": "92233",
                          "eff": separation_efficiency
                        },
                        {
                          "comp": "92234",
                          "eff": separation_efficiency
                        },
                        {
                          "comp": "92235",
                          "eff": separation_efficiency
                        },
                        {
                          "comp": "92236",
                          "eff": separation_efficiency
                        },
                        {
                          "comp": "92238",
                          "eff": separation_efficiency
                        }]
                      }
                    }
                  }
                ]
              }
            }}
          },
          {
            "name": "NuclearWasteSink",
            "config": {"Sink": {
              "in_commods": {"val": ["NuclearWaste"]},
              "max_inv_size": 1e299,
              "capacity": 1e299
            }}
          },
          {
            "name": "PlutoniumSink",
            "config": {"Sink": {
              "in_commods": {"val": ["WeapongradePu"]},
              "max_inv_size": 1e299,
              "capacity": 1e299
            }}
          },
          {
            "name": "WeapongradeUSink",
            "config": {"Sink": {
            "in_commods": {"val": ["EnrichedU"]},
            "recipe_name": "WeapongradeURecipe",
             "max_inv_size": 1e299,
             "capacity": 1e299
            }}
          },
          {
            "name": "DepletedNaturalUSink",
            "config": {"Sink": {
              "in_commods": {"val": ["DepletedNaturalU"]},
              "max_inv_size": 1e299,
              "capacity": 1e299
            }}
          },
          {
            "name": "DepletedReprocessedUSink",
            "config": {"Sink": {
              "in_commods": {"val": ["DepletedReprocessedU"]},
              "max_inv_size": 1e299,
              "capacity": 1e299
            }}
          },
          {
            "name": "SpentUSink",
            "config": {"Sink": {
              "in_commods": {"val": ["SpentU"]},
              "max_inv_size": 1e299,
              "capacity": 1e299
            }}
          }
        ]}

    return d
