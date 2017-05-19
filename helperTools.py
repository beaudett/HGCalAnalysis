#!/usr/bin/env python
import ROOT, math
import numpy as np

def addDataPoint(data,key,point):
    if key in data: data[key].append(point)
    else: data[key] = [point]

def getHisto(values, hname = "hist", htitle = "hist"):

    nbins = 100
    # detect value type (1D/2D)
    #if "tuple" in type(values[0]):
    if isinstance(values[0], tuple):
        #hist_type = "2d"

        # define histo
        xmin = min([val[0] for val in values])
        xmax = max([val[0] for val in values])
        ymin = min([val[1] for val in values])
        ymax = max([val[1] for val in values])
        hist = ROOT.TH2F(hname,htitle,nbins,xmin,xmax,nbins,ymin,ymax)

        # fill
        for val in values: hist.Fill(val[0],val[1])

    else:
        #hist_type = "1d"

        # define histo
        xmin,xmax  = min(values), max(values)
        hist = ROOT.TH1F(hname,htitle,nbins,xmin,xmax)

        # fill
        for val in values: hist.Fill(val)

    #if "int" in type(value[0]) or "float" in type(value[0]):
    #else:
    #    print("Unknown type %s" % type(value[0]))
    ROOT.SetOwnership(hist,0)
    return hist

def compDeltaVar(items, var = "x", prop = "std"):

    values = [getattr(item,var)() for item in items]
    #values = [item.x() for item in items]
    #print values

    res = getattr(np,prop)(values)
    return res
