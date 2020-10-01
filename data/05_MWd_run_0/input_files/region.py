import sys
sys.path.append(
    "/Users/test/Uni/Masterarbeit/thesis_simulations/input/test_transfers/"
)

import institution

def region(burnup):
    inst = institution.institution(burnup)
    d = {"region": [
          {
            "name": "MyRegion",
            "config": {"NullRegion": None},
            "institution": inst["institution"]
          }
        ]}

    return d
