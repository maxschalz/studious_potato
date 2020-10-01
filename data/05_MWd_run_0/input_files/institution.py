import sys
sys.path.append(
    "/Users/test/Uni/Masterarbeit/thesis_simulations/input/test_transfers/"
)

import facility

def institution(burnup):
    fac = facility.facility(burnup)
    d = {"institution": [
           {
             "name": "MyInstitution",
             "config": {"NullInst": None},
             "initialfacilitylist": {
               "entry": [{"number": 1, "prototype": f["name"]} 
                        for f in fac["facility"]]
             }
           }
        ]}
    return d
