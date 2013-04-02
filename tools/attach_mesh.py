#!/usr/bin/env python

import sys
import json

if __name__ == "__main__":

    # open .kem file
    with open(sys.argv[-2], 'rb') as fd:
        meshdata = json.load(fd)

    # open .kep file
    with open(sys.argv[-3], 'rb') as fd:
        data = json.load(fd)

    # attach
    data['mesh'] = meshdata
        
    # dump in .kei file
    with open(sys.argv[-1], 'wb') as fd: 
#        json.dump(data, fd)
        json.dump(data, fd, indent=2, sort_keys=True)

