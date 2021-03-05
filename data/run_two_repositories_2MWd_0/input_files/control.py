def control():
    # Durations in seconds:
    # 2629846  month
    # 86400    day
    # 3600     hour
    d = {"control": {
          "startyear": 2020,
          "startmonth": 1,
          "duration": 730,
          "dt": 86400, 
          "simhandle": "Two years 2MWd and inventory and decay",
          "explicit_inventory": True,
          "explicit_inventory_compact": True,
          "solver": {"config": {"coin-or": {"verbose": False}}},
          "decay": "lazy"
        }}

    return d
