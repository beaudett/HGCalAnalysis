#!/usr/bin/env python
import sys
import ROOT, math
import numpy as np
from NtupleDataFormat import HGCalNtuple
from helperTools import *
# The purpose of this file is to demonstrate mainly the objects
# that are in the HGCalNtuple

max_events = 100

z_half = +1
minE = .5
min_pt = 5

min_fbrem = 0.9
#max_dR = 0.3
max_dR = 0.3

def main(fname = "hgcalNtuple-El15-100_noReClust.root"):
    ntuple = HGCalNtuple(fname)

    foutname = fname.replace(".root","_plots.root")
    tfile = ROOT.TFile(foutname,"recreate")

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

    hist_data = {}

    # store rechits
    good_cluster_rechits = {}
    bad_cluster_rechits = {}

    for event in ntuple:
        #if stop_run: break
        if tot_nevents >= max_events: break

        if tot_nevents % 100 == 0: print("Event %i" % tot_nevents)
        # print "Event", event.entry()
        tot_nevents += 1

        genParts = event.genParticles()
        #tot_genpart += len(genParts)

        found_part = False
        for part in genParts:

            if part.eta() * z_half < 0: continue
            if part.gen() < 1: continue
            if not part.reachedEE(): continue
            if part.fbrem() < min_fbrem: continue

            found_part = True
            tot_genpart += 1
            break

        if not found_part: continue

        multiClusters = event.multiClusters()
        #tot_multiclus += len(multiClusters)

        # 2d clusters in layers
        layerClusters = event.layerClusters()
        tot_cluster2d += len(layerClusters)

        # rechits
        recHits = event.recHits()
        tot_rechit += len(recHits)

        for i_mcl, multicl in enumerate(multiClusters):
            if multicl.energy() < minE: continue
            if multicl.pt() < min_pt: continue
            if multicl.eta() * part.eta() < 0: continue
            if len(multicl.cluster2d()) < 3: continue

            dR =  math.hypot(multicl.eta()-part.eta(), multicl.phi()-part.phi())

            if dR > max_dR: continue
            tot_multiclus += 1

            addDataPoint(hist_data,"part_mcl_dR",dR)
            #addDataPoint(hist_data,"mcl_axisZ",multicl.pcaAxisZ())

            ## play with cluster layers
            clusters = [layerClusters[clust_idx] for clust_idx in multicl.cluster2d()]

            #gr_rh_XYZ = ROOT.TGraph2D(); gr_rh_XYZ.SetMarkerStyle(20)

            ## Rechit plots
            rechits = []
            for cluster in clusters:
                rechits += [recHits[rh_idx] for rh_idx in cluster.rechits()]

            rh_coords = [(rh.x(),rh.y(),rh.z()) for rh in rechits]

            gr_rh_XYZ = ROOT.TGraph2D(); gr_rh_XYZ.SetMarkerStyle(20)
            gr_rh_XZ = ROOT.TGraph(); gr_rh_XZ.SetMarkerStyle(20)
            gr_rh_YZ = ROOT.TGraph(); gr_rh_YZ.SetMarkerStyle(20)

            hXZ = ROOT.TH2F("hXZ_event_%i"%event.event(),"XZ rechits",120,320,350,150,-150,150)

            #for i,rh in enumerate(clusters):
            for i,rh in enumerate(rechits):

                if rh.layer() > 28: continue
                #if rh.detid() > 1.2 * 10e9: continue

                gr_rh_XYZ.SetPoint(i,rh.z(),rh.x(),rh.y())
                gr_rh_XZ.SetPoint(i,rh.z(),rh.x())
                gr_rh_YZ.SetPoint(i,rh.z(),rh.y())

                hXZ.Fill(rh.z(),rh.x(),rh.energy())

                if abs(multicl.pcaAxisZ()) > 0.5:
                    dr = math.hypot(multicl.eta()-rh.eta(), multicl.phi()-rh.phi())
                    addDataPoint(hist_data,"gd_mcl_rh_dR",dr)
                    drho = math.hypot(multicl.slopeX()-rh.x(), multicl.slopeY()-rh.y())
                    addDataPoint(hist_data,"gd_mcl_rh_drho",drho)
                    #addDataPoint(hist_data,"gd_mcl_rh_drho_E",(drho,rh.layer()))
                    addDataPoint(hist_data,"gd_mcl_rh_drho_E",(drho,rh.z()))
                else:
                    dr = math.hypot(multicl.eta()-rh.eta(), multicl.phi()-rh.phi())
                    addDataPoint(hist_data,"bad_mcl_rh_dR",dr)
                    drho = math.hypot(multicl.slopeX()-rh.x(), multicl.slopeY()-rh.y())
                    addDataPoint(hist_data,"bad_mcl_rh_drho",drho)
                    #addDataPoint(hist_data,"bad_mcl_rh_drho_E",(drho,rh.layer()))
                    addDataPoint(hist_data,"bad_mcl_rh_drho_E",(drho,rh.z()))

            grtitle = "axisZ %0.2f, event %i" % (abs(multicl.pcaAxisZ()),event.event())

            if abs(multicl.pcaAxisZ()) > 0.5:
                #good_cluster_rechits.append(rh_coords)
                good_cluster_rechits[event.event()] = rh_coords

                #grtitle = "Good axisZ %f, event %i" % (abs(multicl.pcaAxisZ())event.event())
                grtitle = "Good " + grtitle
                grname = "grxyz_good_event_%i" % event.event()

            else:
                bad_cluster_rechits[event.event()] = rh_coords
                #bad_cluster_rechits.append(rh_coords)

                #grtitle = "Bad axisZ, event %i" % event.event()
                grtitle = "Bad " + grtitle
                grname = "grxyz_bad_event_%i" % event.event()

            gr_rh_XYZ.SetTitle(grtitle)
            gr_rh_XYZ.SetName(grname)

            cname = grname.replace("gr","c")
            canv = ROOT.TCanvas(cname,"Event",1000,800)
            canv.Divide(2,2)

            canv.cd(1); gr_rh_XYZ.Draw("p")
            canv.cd(2); gr_rh_XZ.Draw("ap")
            canv.cd(3); gr_rh_YZ.Draw("ap")
            canv.cd(4); hXZ.Draw("colz")

            canv.Update()
            tfile.cd()
            canv.Write()
            #q = raw_input("Cont..")
        #break

    tfile.Close()

    '''
    print "Good clusters"
    for event,rechits in good_cluster_rechits.iteritems():
        print len(rechits)

    print "Bad clusters"
    for event,rechits in bad_cluster_rechits.iteritems():
        print len(rechits)
    '''
    hists = []
    for data_name in hist_data:
        print("Plotting hist for data: %s" %data_name)
        canv = ROOT.TCanvas("canv_" + data_name,data_name,800,600)
        hist = getHisto(hist_data[data_name],data_name,data_name)
        hist.Draw("colz")
        canv.Update()
        #hists.append(hist)
        ROOT.SetOwnership(canv,0)
    q = raw_input("exit")

if __name__ == "__main__":

    #if '-b' in sys.argv: sys.argv = [sys.argv[0]]

    if len(sys.argv) > 1:
        if '-b' in sys.argv: main()
        else:
            fname = sys.argv[1]
            print '# Input file is', fname
            main(fname)
    else:
        print("No input files given!")
        main()
