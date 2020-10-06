import sys
sys.path.append(
    "/Users/test/Uni/Masterarbeit/thesis_simulations/input/test_transfers/"
)

def archetypes():
    d = {"archetypes": {
           "spec": [
             {"lib": "agents", "name": "NullInst"},
             {"lib": "agents", "name": "NullRegion"},
             {"lib": "cycamore", "name": "Reactor"},
             {"lib": "cycamore", "name": "Separations"},
             {"lib": "cycamore", "name": "Sink"},
             {"lib": "cycamore", "name": "Source"},
             {"lib": "cycamore", "name": "Storage"},
             {"lib": "misoenrichment", "name":"MIsoEnrich"}
           ]
         }}

    return d
