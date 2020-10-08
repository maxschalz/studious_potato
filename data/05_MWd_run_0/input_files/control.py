import sys
sys.path.append(
    "/Users/test/Uni/Masterarbeit/thesis_simulations/input/test_transfers/"
)

def control():
    # Durations in seconds:
    # 2629846  month
    # 86400    day
    # 3600     hour
    d = {"control": {
          "startyear": 2020,
          "startmonth": 1,
          "duration": 365,
          "dt": 86400, 
          "simhandle": "Test flows",
          "explicit_inventory": True,
          "explicit_inventory_compact": True,
          "solver": {"config": {"coin-or": {"verbose": False}}},
          "decay": "lazy"
        }}

    return d
