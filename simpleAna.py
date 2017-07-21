#!/usr/bin/env python
import sys
import ROOT, math, os
import numpy as np
from NtupleDataFormat import HGCalNtuple
import collections
from helperTools import *

# The purpose of this file is to demonstrate mainly the objects
# that are in the HGCalNtuple

max_events = 1000

minE = 0
min_pt = 0

min_fbrem = -1
#max_dR = 0.3
max_dR = 0.3
twopi = 2*math.pi
histDict = {}
histToShow = {}

def addToDict(histo,order=0):
    histoName = histo.GetName()
    histDict[histoName]=histo
    if order > 0:
        if not order in histToShow:        
            histToShow[order]=[]
        histToShow[order].append(histoName)

def anaTree(opts):
    #ntuple = HGCalNtuple("/Users/clange/CERNBox/partGun_PDGid211_x120_E80.0To80.0_NTUP_9.root")
    ntuple = HGCalNtuple(opts.infile)
#    ntuple = HGCalNtuple("hgcalNtuple-pca.root")
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

#    addToDict(ROOT.TH1F("h_mclust_dR","multi dR; dR(p,mclust)",100,0,3.2))
    addToDict(ROOT.TH1F("h_clust_dPhi","h_clust_dPhi;dPhi(mclust,clust)",100,-0.1,0.1),1)    
    addToDict(ROOT.TH1F("h_clust_dEta","h_clust_dEta;dEta(mclust,clust)",100,-0.1,0.1),1)
    addToDict(ROOT.TH1F("h_clust_dr","h_clust_dr",80,0,4),1)
    addToDict(ROOT.TH1F("h_clust_drtMin","h_clust_drtMin",100,-2,2),1)
    addToDict(ROOT.TH1F("h_clust_dPhirtMin","h_clust_dPhirtMin",100,-0.1,0.1),1)
    addToDict(ROOT.TH2F("h_E_vs_gen","h_E_vs_gen",100,0,15,100,0,15),0)
    addToDict(ROOT.TH1F("h_clust_dr_newdist","association distance",100,0,10),1)
    
#    addToDict(ROOT.TH1F("h_clust_dR","h_clust_dR;dR(mclust,clust)",100,0,0.2))
    addToDict(ROOT.TH1F("h_clust_mult_tot_p","Pt tot, eta>0 ",100,0,50),0)
    addToDict(ROOT.TH1F("h_clust_mult_tot_n","Pt tot, eta<0 ",100,0,50),0)
    addToDict(ROOT.TH1F("h_clust_mult_siguu","#sigma_{uu} ",75,0,15),4)
    addToDict(ROOT.TH1F("h_clust_mult_sigvv","#sigma_{vv} ",75,0,15),4)
    addToDict(ROOT.TH1F("h_clust_mult_siguu_lowBrem","#sigma_{uu} - low Brem; #sigma_{uu}",75,0,15),6)
    addToDict(ROOT.TH1F("h_clust_mult_sigvv_lowBrem","#sigma_{vv} - low Brem; #sigma_{vv}",75,0,15),7)
    addToDict(ROOT.TH1F("h_clust_mult_siguu_highBrem","#sigma_{uu} - high Brem; #sigma_{uu}",75,0,15),6)
    addToDict(ROOT.TH1F("h_clust_mult_sigvv_highBrem","#sigma_{vv} -high Brem; #sigma_{vv}",75,0,15),7)
    addToDict(ROOT.TProfile("h_clust_mult_rad_siguu","#sigma_{uu} vs r;r;#sigma_{uu}",10,0.5,10.5),5)
    addToDict(ROOT.TProfile("h_clust_mult_rad_sigvv","#sigma_{vv} vs r;r;#sigma_{vv}",10,0.5,10.5),5)
              
#    addToDict(ROOT.TH2F("h_mclust_dR_dZ","multi dR; dR(p,mclust);dZ (clusts)",100,0,3.2,100,0,20))
 
    print str(len(histDict))+" Histos created"

    for event in ntuple:
        if tot_nevents >= max_events: break
        tot_multiclus_pos = 0
        tot_multiclus_neg = 0
        multiClusters = event.multiClusters()
        for  i_mcl,multicl in enumerate(multiClusters):
            if ( multicl.eta() > 0.):
                tot_multiclus_pos += multicl.pt()
            elif (multicl.eta() < 0.):
                tot_multiclus_neg += multicl.pt()

        if tot_nevents % 100 == 0: print("Event %i" % tot_nevents)
        # print "Event", event.entry()
        tot_nevents += 1

        genParts = event.genParticles()
        # 2d clusters in layers
        layerClusters = event.layerClusters()
        tot_cluster2d += len(layerClusters)

        # rechits
        recHits = event.recHits()
        tot_rechit += len(recHits)

        #tot_genpart += len(genParts)
#        print "Part"+str(len(genParts))
        found_part = False
        for part in genParts:
            if not abs(part.pid()) == 11: continue
            if part.gen() < 1: continue
            if not part.reachedEE()==2 : continue
            
            found_part = True
            
            tot_genpart += 1

            if not found_part: continue
                        
            dRMin = 999.
            drMin = 9999.
            idrMin = -1
            dPhidrMin = 999.
            dEtadrMin = 999.
            layer=10
            drtMin = 999.
            irtMin = -1
            dPhirtMin = 999.
            for i_mcl, multicl in enumerate(multiClusters):
                #            print str(multicl.energy())+" "+ str(multicl.pt())+" "+str(len(multicl.cluster2d()))
                if multicl.energy() < minE: continue
                if multicl.pt() < min_pt: continue
                if multicl.eta() * part.eta() < 0: continue
                

                if len(multicl.cluster2d()) < 3: continue
                
                #if abs(multicl.pcaAxisZ()) < 0.4: continue
                dphi = multicl.phi()-part.phi()
                deta = multicl.eta()-part.eta()
                drt = math.hypot(part.posx()[layer],part.posy()[layer])-math.hypot(multicl.pcaPosX(),multicl.pcaPosY())
                if math.fabs(drt) < math.fabs(drtMin):
                    drtMin = drt
                    irtMin = i_mcl
                    dPhirtMin = math.fmod(math.atan2(multicl.pcaPosY(),multicl.pcaPosX())-math.atan2(part.posy()[layer],part.posx()[layer]),twopi)
                                            
                if (dphi > math.pi):
                    dphi -= twopi
                elif (dphi < -math.pi):
                    dphi += twopi
                    
                dR =  math.hypot(multicl.eta()-part.eta(),dphi )
                if dR < dRMin:
                    dRMin = dR

                    dr = (part.posx()[layer]-multicl.pcaPosX())* (part.posx()[layer]-multicl.pcaPosX()) + (part.posy()[layer]-multicl.pcaPosY())*(part.posy()[layer]-multicl.pcaPosY())
                    if dr<drMin:
                        drMin = dr 
                        idrMin = i_mcl
                        #                dPhidrMin = dphi
                        dPhidrMin = math.atan2(multicl.pcaPosY(),multicl.pcaPosX())-math.atan2(part.posy()[layer],part.posx()[layer])
                        if (dPhidrMin > math.pi):
                            dPhidrMin -= twopi
                        elif (dPhidrMin < -math.pi):
                            dPhidrMin += twopi
                            dEtadrMin = deta
            drMin = math.sqrt(drMin)            
            histDict['h_clust_dPhi'].Fill(dPhidrMin)
            histDict['h_clust_dEta'].Fill(dEtadrMin)
            histDict['h_clust_dr'].Fill(drMin)
            histDict['h_clust_drtMin'].Fill(drtMin)
            if math.fabs(drtMin) < 0.6:
                histDict['h_clust_dPhirtMin'].Fill(dPhirtMin)
                if math.fabs(dPhirtMin) < 0.01:
                    pout = (1.-part.fbrem())*part.pt()
                    histDict['h_E_vs_gen'].Fill(pout,multiClusters[irtMin].pt())
                    dist = math.sqrt((part.posx()[layer]-multiClusters[irtMin].pcaPosX())* (part.posx()[layer]-multiClusters[irtMin].pcaPosX()) + (part.posy()[layer]-multiClusters[irtMin].pcaPosY())*(part.posy()[layer]-multiClusters[irtMin].pcaPosY()))
                    histDict['h_clust_dr_newdist'].Fill(dist)

            histDict['h_clust_mult_tot_p'].Fill(tot_multiclus_pos)
            histDict['h_clust_mult_tot_n'].Fill(tot_multiclus_neg)
            
            if(drMin < 2 and multiClusters[idrMin].NLay() > 2 ):
                histDict['h_clust_mult_siguu'].Fill(multiClusters[idrMin].siguu())
                histDict['h_clust_mult_sigvv'].Fill(multiClusters[idrMin].sigvv())
                    
                if part.fbrem() < 0.5:
                    histDict['h_clust_mult_siguu_lowBrem'].Fill(multiClusters[idrMin].siguu())
                    histDict['h_clust_mult_sigvv_lowBrem'].Fill(multiClusters[idrMin].sigvv())
                else:
                    histDict['h_clust_mult_siguu_highBrem'].Fill(multiClusters[idrMin].siguu())
                    histDict['h_clust_mult_sigvv_highBrem'].Fill(multiClusters[idrMin].sigvv())
                    
                if ("radsig" not in options.skipflags):
                    for radius,siguu,sigvv in zip(multiClusters[idrMin].sigrad(),multiClusters[idrMin].radsiguu(),multiClusters[idrMin].radsigvv()):
                        histDict['h_clust_mult_rad_siguu'].Fill(radius,siguu)
                        histDict['h_clust_mult_rad_sigvv'].Fill(radius,sigvv)

    
                    
    print "Processed %d events" % tot_nevents
    flushHistos(opts.outfile,histDict)
    if not options.plotdir == "":
        makePdfs(options.plotdir,histDict)

    if opts.display:
        displayPlots(histDict)

#    q = raw_input("exit")

def flushHistos(histoFile,histoDict):    
    output = ROOT.TFile(histoFile,"RECREATE")
    for name,histo in histDict.items():
        histo.Write()

    output.Write()
    output.Close()

def displayPlots(histoDict):
    print "Display plots"

    od = collections.OrderedDict(sorted(histToShow.items()))
    for key,names in od.iteritems():
        for name in names:
            c1 = ROOT.TCanvas(name,name)
            histoDict[name].Draw()
            c1.Draw()
            c1.Update()

def makePdfs(directory,histoDict):
    if not os.path.exists(options.plotdir): os.makedirs(options.plotdir)
    for name,histo in histoDict.items():
        c1 = ROOT.TCanvas(False)
        c1.SetBatch(False)
        histo.Draw()
        c1.Update()
        c1.Print(options.plotdir+"/"+name+".pdf")

if __name__ == "__main__":
    from optparse import OptionParser
    parser = OptionParser()
    
    parser.usage = '%prog [options]'
    parser.description="""
    Analyse HGCAL flat tuple
    """
    
    # General options
#    parser.add_option("-b","--batch", dest="batch",default=False, action="store_true", help="Batch mode")
    parser.add_option("-n","--maxEve","--maxEvents", "--maxEntries",dest="maxEntries", default=-1,  type="int",    help="Maximum entries to analyze")
    parser.add_option("-v","--verbose",  dest="verbose",  default=1,  type="int",    help="Verbosity level (0 = quiet, 1 = verbose, 2+ = more)")

    # File options
    parser.add_option("-f","--file", dest="infile",default="hgcalNtuple-pca.root",help="Input File name(s)")
#    parser.add_option("-t","--tname", dest="treename",default="ana/hgc",help="Tree name")

    # Plot options
    parser.add_option("-d","--display", dest="display",default=False,help="Display plots",action="store_true")
    parser.add_option("-o","--outfile", dest="outfile",default="plots.root",help="Plot file")
    parser.add_option("-p","--pdir","--plotdir", dest="plotdir",default="",help="Plot directory")
    parser.add_option("--skipflag",dest="skipflags",action="append",help="Run-time flags",default=[''])
#    parser.add_option("-l","--liTerm2",dest="listplots",help="Show plots in iTerm2",default=False,action="store_true")

    # Read options and args
    (options,args) = parser.parse_args()

    # evaluate options


    ###############################################################
    ### Start running
    ###############################################################

    print("Starting")
    anaTree(options)
