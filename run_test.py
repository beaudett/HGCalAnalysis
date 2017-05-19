#!/usr/bin/env python
import ROOT, math
import numpy as np
from NtupleDataFormat import HGCalNtuple

# The purpose of this file is to demonstrate mainly the objects
# that are in the HGCalNtuple

z_half = +1
minE = .5
min_fbrem = 0.6
#max_dR = 0.3
max_dR = 3

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
    return hist

def compDeltaVar(items, var = "x", prop = "std"):

    values = [getattr(item,var)() for item in items]
    #values = [item.x() for item in items]
    #print values

    res = getattr(np,prop)(values)
    return res

def main():
    #ntuple = HGCalNtuple("/Users/clange/CERNBox/partGun_PDGid211_x120_E80.0To80.0_NTUP_9.root")
    ntuple = HGCalNtuple("hgcalNtuple-pca-1000.root")
    #ntuple = HGCalNtuple("../hgcalNtuple-pca.root")

    tot_nevents = 0
    tot_genpart = 0
    tot_rechit = 0
    tot_rechit_raw = 0
    tot_cluster2d = 0
    tot_multiclus = 0
    tot_simcluster = 0
    tot_pfcluster = 0
    tot_calopart = 0
    tot_track = 0

    h_mclust_dR = ROOT.TH1F("h_mclust_dR","multi dR; dR(p,mclust)",100,0,3.2)

    h_clust_dEta = ROOT.TH1F("h_clust_dEta","h_clust_dEta;dEta(mclust,clust)",100,-0.1,0.1)
    h_clust_dPhi = ROOT.TH1F("h_clust_dPhi","h_clust_dPhi;dPhi(mclust,clust)",100,-0.1,0.1)
    h_clust_dR = ROOT.TH1F("h_clust_dR","h_clust_dR;dR(mclust,clust)",100,0,0.2)

    h_mclust_dR_dZ = ROOT.TH2F("h_mclust_dR_dZ","multi dR; dR(p,mclust);dZ (clusts)",100,0,3.2,100,0,20)
    #h_mclust_dXdY = ROOT.TH2F("h_mclust_dXdY","dX/dY; layer; dX/dY",28,0,28,100,0,30)
    h_mclust_dX = ROOT.TH2F("h_mclust_dX","dX; layer; dX",28,1,29,100,0,10)
    h_mclust_dY = ROOT.TH2F("h_mclust_dY","dY; layer; dY",28,1,29,100,0,10)

    ## data storages for hists
    hist_data = {}

    for event in ntuple:
        if tot_nevents >= 100: break

        if tot_nevents % 100 == 0: print("Event %i" % tot_nevents)
        # print "Event", event.entry()
        tot_nevents += 1
        '''
        if (ntuple.hasRawRecHits()):
            recHitsRaw = event.recHits("rechit_raw")
            tot_rechit_raw += len(recHitsRaw)
        layerClusters = event.layerClusters()
        tot_cluster2d += len(layerClusters)
        simClusters = event.simClusters()
        tot_simcluster += len(simClusters)
        pfClusters = event.pfClusters()
        tot_pfcluster += len(pfClusters)
        pfClusters = event.pfClusters()
        tot_pfcluster += len(pfClusters)
        caloParts = event.caloParticles()
        tot_calopart += len(caloParts)
        tracks = event.tracks()
        tot_track += len(tracks)
        '''
        genParts = event.genParticles()
        #tot_genpart += len(genParts)

        has_brem = False
        found_part = False
        for part in genParts:

            if part.eta() * z_half < 0: continue
            if part.gen() < 1: continue
            if not part.reachedEE(): continue
            if part.fbrem() < min_fbrem: continue

            found_part = True
            tot_genpart += 1

            '''
            #print part.fbrem()
            has_brem = True
            tot_genpart += 1

            print part.energy(), part.pt(), part.gen(), part.fbrem(), part.eta(), (part.eta() * z_half < 0)
            # store particle in part
            '''
        #if has_brem: continue
        if not found_part: continue

        multiClusters = event.multiClusters()
        #tot_multiclus += len(multiClusters)

        # 2d clusters in layers
        layerClusters = event.layerClusters()
        tot_cluster2d += len(layerClusters)

        # rechits
        recHits = event.recHits()
        tot_rechit += len(recHits)


        for multicl in multiClusters:
            if multicl.energy() < minE: continue
            if multicl.eta() * part.eta() < 0: continue

            if len(multicl.cluster2d()) < 3: continue

            dR =  math.hypot(multicl.eta()-part.eta(), multicl.phi()-part.phi())
            h_mclust_dR.Fill(dR)

            tot_multiclus += 1

            if dR > max_dR: continue

            addDataPoint(hist_data,"clust_dR",dR)
            '''
            print "Here", multicl.energy()
            print multicl.eigenVal1(), multicl.eigenVal2(), multicl.eigenVal3()
            #if multicl.pcaAxisX() < multicl.pcaAxisY():
            print multicl.pcaAxisX(), multicl.pcaAxisY(), multicl.pcaAxisZ()
            print multicl.pcaPosX(), multicl.pcaPosY(), multicl.pcaPosZ()

            print multicl.phi(), multicl.eta()
            '''

            '''
            for clust_idx in multicl.cluster2d():
                clust = layerClusters[clust_idx]
                if clust.layer() > 28: continue

                #print clust.layer(),
                #print clust, layerClusters[clust].phi(), layerClusters[clust].eta()
                clust_dR =  math.hypot(multicl.eta()-clust.eta(), multicl.phi()-clust.phi())
                h_clust_dR.Fill(clust_dR)
                h_clust_dEta.Fill(clust.phi()-multicl.phi())
                h_clust_dPhi.Fill(clust.eta()-multicl.eta())
            '''

            ## play with cluster layers
            clusters = [layerClusters[clust_idx] for clust_idx in multicl.cluster2d()]
            #print len(clusters), " clusters in multi"
            #print compDeltaVar(clusters,'z',"mean"),
            dZ = compDeltaVar(clusters,'z',"std")
            #print dZ
            h_mclust_dR_dZ.Fill(dR,dZ)

            #continue
            #if len(clusters) < 25: continue

            cluster_layers = [[] for layer in range(28)]
            rechits_layers = [[] for layer in range(28)]
            #cluster_layers = [28*[]]
            #rechits_layers = [28*[]]

            for cluster in clusters:
                if cluster.layer() > 28: continue
                #if cluster.energy() < minE/100.: continue
                cluster_layers[cluster.layer()-1].append(cluster)
                # rechits in layer
                rechits_layers[cluster.layer()-1]+=[recHits[rh_idx] for rh_idx in cluster.rechits()]

            for lay,rhits in enumerate(rechits_layers):
                if len(rhits) == 0: continue
                dX = compDeltaVar(rhits,'x',"std")
                dY = compDeltaVar(rhits,'y',"std")

                h_mclust_dX.Fill(lay+1,dX)
                h_mclust_dY.Fill(lay+1,dY)

            '''
            for lay,clusts in enumerate(cluster_layers):
                if len(clusts) == 0: continue
                dX = compDeltaVar(clusts,'x',"std")
                dY = compDeltaVar(clusts,'y',"std")

                h_mclust_dX.Fill(lay,dX)
                h_mclust_dY.Fill(lay,dY)


                print("%i - %i , dX = %f" %(lay, len(clusts), dX))

                #for prop in ["std","mean"]:
                #    print prop
                #    for var in ["x","y","z"]:
                #        print var, compDeltaVar(clusts, var, prop),

            #print
            '''


        # for genPart in genParts:
        #     print tot_nevents, "genPart pt:", genPart.pt()

    print "Processed %d events" % tot_nevents
    print "On average %f generator particles" % (float(tot_genpart) / tot_nevents)
    print "On average %f reconstructed hits" % (float(tot_rechit) / tot_nevents)
    print "On average %f raw reconstructed hits" % (float(tot_rechit_raw) / tot_nevents)
    print "On average %f layer clusters" % (float(tot_cluster2d) / tot_nevents)
    print "On average %f multi-clusters" % (float(tot_multiclus) / tot_nevents)
    print "On average %f sim-clusters" % (float(tot_simcluster) / tot_nevents)
    print "On average %f PF clusters" % (float(tot_pfcluster) / tot_nevents)
    print "On average %f calo particles" % (float(tot_calopart) / tot_nevents)
    print "On average %f tracks" % (float(tot_track) / tot_nevents)

    canv_dR = ROOT.TCanvas("canv_dR","canv",800,600)
    h_mclust_dR_dZ.Draw("colz")
    canv_dR.Update()

    canv_dXY = ROOT.TCanvas("canv_dXY","canv",1000,600)
    canv_dXY.Divide(2,1)

    canv_dXY.cd(1)
    h_mclust_dX.Draw("colz")
    canv_dXY.cd(2)
    h_mclust_dY.Draw("colz")

    canv_dXY.Update()

    '''
    canv.SetLogy()

    canv.cd(1)
    h_mclust_dR.Draw()
    #h_clust_dEta.Draw()
    canv.cd(2)
    #h_clust_dPhi.Draw()
    h_clust_dR.Draw()
    '''

    #print hist_data
    for data_name in hist_data:
        canv = ROOT.TCanvas("canv_" + data_name,data_name,800,600)
        hist = getHisto(hist_data[data_name])
        hist.Draw()
        canv.Update()
    q = raw_input("")

if __name__ == "__main__":
    main()
