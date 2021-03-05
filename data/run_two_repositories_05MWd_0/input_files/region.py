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
