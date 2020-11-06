import sys
sys.path.append(
    "/Users/test/Uni/Masterarbeit/thesis_simulations/input/test_transfers/"
)

def commodity():
    d = {"commodity": [
          {"name": "WeapongradePu", "solution_priority": 5},
          {"name": "SpentFuel", "solution_priority": 5},
          {"name": "SpentReprocessedFuel", "solution_priority": 5},
          {"name": "EnrichedU", "solution_priority": 1}
        ]}

    return d
