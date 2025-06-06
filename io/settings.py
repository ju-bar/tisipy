import json

def load(sfile, format='json', debug=False):
    with open(sfile, 'r') as f:
        d = json.load(f)
        if debug:
            if "title" in d:
                print("loaded " + d["title"] + " from file " + sfile)
        return d
    return {}