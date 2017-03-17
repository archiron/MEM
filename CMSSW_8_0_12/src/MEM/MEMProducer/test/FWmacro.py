from sys import argv
argv.append( '-b-' )
import ROOT
ROOT.gROOT.SetBatch(True)
argv.remove( '-b-' )

ROOT.gSystem.Load("libFWCoreFWLite.so");
ROOT.gSystem.Load("libDataFormatsFWLite.so");
ROOT.FWLiteEnabler.enable()

from ROOT import *
from math import log10

from DataFormats.FWLite import Handle, Events
from optparse import Values
options = Values()
options.inputFiles = [ argv[1] ]
options.secondaryInputFiles = argv[2:] 
options.maxEvents = -1
events = Events([],options=options)

def createPicture(histo1, histo2, fileName):

    print fileName
    c1.Clear()
    c1.SetLogy(0)
    histo2.Draw()
    c1.Update()
    gMax2 = ROOT.gPad.GetUymax()
    print "    histo 2 max : ", ROOT.gPad.GetUymax(), gMax2

    c1.Clear()
    histo1.Draw()
    c1.Update()
    gMax1 = ROOT.gPad.GetUymax()
    print "    histo 1 max : ", ROOT.gPad.GetUymax(), gMax1

    var_1 = log10( abs(gMax1 - gMax2) )
    print "log = %8.2f" % var_1
    if (var_1 > 1.5 ): # and histo1.GetMaximum() > 0
        c1.SetLogy(1)
    
    c1.Clear()
    if (gMax2 >= gMax1):
        print "    gMax2 > gMax1"
        histo2.SetLineColor(4)       #red
        histo2.SetMinimum(0.5)
        histo2.Draw()
        histo1.SetLineColor(2)       #blue
        histo1.SetMinimum(0.5)
        histo1.Draw("histsames") # histsames
        c1.Update()
        statBox2 = histo2.GetListOfFunctions().FindObject("stats")
        statBox2.SetTextColor(kBlue)
        y1 = statBox2.GetY1NDC()
        y2 = statBox2.GetY2NDC()
        statBox2.SetY1NDC(2*y1-y2) ;
        statBox2.SetY2NDC(y1) ;
        c1.Update()
        statBox1 = histo1.GetListOfFunctions().FindObject("stats")
        statBox1.SetTextColor(kRed)
        c1.Update()
    else:
        print "    gMax1 > gMax2"
        histo1.SetLineColor(2)       #blue
        histo1.SetMinimum(0.5)
        histo1.Draw() # 
        histo2.SetLineColor(4)       #red
        histo2.SetMinimum(0.5)
        histo2.Draw("histsames")
        c1.Update()
        statBox2 = histo2.GetListOfFunctions().FindObject("stats")
        statBox2.SetTextColor(kBlue)
        y1 = statBox2.GetY1NDC()
        y2 = statBox2.GetY2NDC()
        statBox2.SetY1NDC(2*y1-y2) ;
        statBox2.SetY2NDC(y1) ;
        c1.Update()
        statBox1 = histo1.GetListOfFunctions().FindObject("stats")
        statBox1.SetTextColor(kRed)
        c1.Update()
    
    Hist_name = dossier + "/" + fileName
    c1.Print (Hist_name)

def createPicture_2(histo1, histo2, fileName):

    print fileName
    c1.Clear()
    c1.SetLogy(0)
    histo2.Draw()
    histo2.SetStats(1)
    ROOT.gPad.Update()
    gMax2 = ROOT.gPad.GetUymax()
    print "    histo 2 max : ", ROOT.gPad.GetUymax(), gMax2
    statBox2 = histo2.GetListOfFunctions().FindObject("stats")
    histo2.SetLineColor(kBlue)       #blue
    histo2.SetMinimum(0.5)
    histo2.SetMarkerColor(kBlue) 
    statBox2.SetTextColor(kBlue)

    ROOT.gPad.Update()
    histo1.Draw()
    histo1.SetStats(1)
    c1.Update()
    gMax1 = ROOT.gPad.GetUymax()
    print "    histo 1 max : ", ROOT.gPad.GetUymax(), gMax1
    statBox1 = histo1.GetListOfFunctions().FindObject("stats")
    histo1.SetLineColor(kRed)       #red
    histo1.SetMinimum(0.5)
    histo1.SetMarkerColor(kRed) ;
    statBox1.SetTextColor(kRed)

    var_1 = log10( abs(gMax1 - gMax2) )
    print "log = %8.2f" % var_1
    if (var_1 > 1.5 ): # and histo1.GetMaximum() > 0
        c1.SetLogy(1)

    if histo2.GetMaximum() > histo1.GetMaximum():
        histo1.SetMaximum(histo2.GetMaximum()*1.1)
    
    y1 = statBox2.GetY1NDC()
    y2 = statBox2.GetY2NDC()
    statBox2.SetY1NDC(2*y1-y2) 
    statBox2.SetY2NDC(y1) 
    
    histo1.Draw()
    histo2.Draw("histsames")
    c1.Update()

    Hist_File_Name = dossier + "/" + fileName # dossier + 
    c1.Print (Hist_File_Name)

dossier = options.inputFiles[0].split('/')[0]

print "In total there are %d events" % events.size()

# Create histograms, etc.
ROOT.gROOT.SetStyle('Plain') # white background
c1 = ROOT.TCanvas()

muonPtHist_1 = ROOT.TH1F ("muonPt", "muon Pt", 50, 0, 150)
muonPtHist_2 = ROOT.TH1F ("muonPt", "muon Pt", 50, 0, 150)
muonChargeHist_1 = ROOT.TH1F ("muonCharge", "muon Charge", 50, -10., 10.)
muonChargeHist_2 = ROOT.TH1F ("muonCharge", "muon Charge", 50, -10., 10.)
metPtHist_1 = ROOT.TH1F ("metPt", "met Pt", 50, 0, 150)
metPtHist_2 = ROOT.TH1F ("metPt", "met Pt", 50, 0, 150)
metChargeHist_1 = ROOT.TH1F ("metCharge", "met Charge", 50, -10., 10.)
metChargeHist_2 = ROOT.TH1F ("metCharge", "met Charge", 50, -10., 10.)
tauPtHist_1 = ROOT.TH1F ("tauPt", "tau Pt", 50, 0, 150)
tauPtHist_2 = ROOT.TH1F ("tauPt", "tau Pt", 50, 0, 150)
tauChargeHist_1 = ROOT.TH1F ("tauCharge", "tau Charge", 50, -10., 10.)
tauChargeHist_2 = ROOT.TH1F ("tauCharge", "tau Charge", 50, -10., 10.)
jetPtHist_1 = ROOT.TH1F ("jetPt", "jet Pt", 50, 0, 150)
jetPtHist_2 = ROOT.TH1F ("jetPt", "jet Pt", 50, 0, 150)
jetChargeHist_1 = ROOT.TH1F ("jetCharge", "jet Charge", 50, -10., 10.)
jetChargeHist_2 = ROOT.TH1F ("jetCharge", "jet Charge", 50, -10., 10.)
elecPtHist_1 = ROOT.TH1F ("Pt", "electron Pt", 50, 0, 150)
elecPtHist_2 = ROOT.TH1F ("Pt", "electron Pt", 50, 0, 150)
elecChargeHist_1 = ROOT.TH1F ("Charge", "electron Charge", 50, -10., 10.)
elecChargeHist_2 = ROOT.TH1F ("Charge", "electron Charge", 50, -10., 10.)
chiPtHist_2 = ROOT.TH1F ("metPt", "met Pt", 50, 0, 150)
chiChargeHist_2 = ROOT.TH1F ("metCharge", "met Charge", 50, -10., 10.)

muonH, muonN = Handle("std::vector<pat::Muon>"), "slimmedMuons"
metH, metN = Handle("std::vector<pat::MET>"), "slimmedMETs"
tauH, tauN = Handle("std::vector<pat::Tau>"), "slimmedTaus"
jetH, jetN = Handle("std::vector<pat::Jet>"), "slimmedJets"
#eleH, eleN = Handle("std::vector<pat::Electron>"), "slimmedElectrons"

print "boucle enumerate 1"
iTauCount_1 = 0
iMuonCount_1 = 0
for i,event in enumerate(events):
    print "Event : ", i
    event.getByLabel(muonN, muonH)
    event.getByLabel(metN, metH)
    event.getByLabel(tauN, tauH)
    event.getByLabel(jetN, jetH)
#    event.getByLabel(eleN, eleH)
        
    for evt in muonH.product():
        print "Muon      pT %8.2f " % evt.pt()
        iMuonCount_1 += 1
        muonPtHist_1.Fill( evt.pt() )
        muonChargeHist_1.Fill( evt.charge() )

    for evt in metH.product():
        print "MET       pT %8.2f " % evt.pt()
        print "MET uncor pT %8.2f " % evt.uncorPt()
        metPtHist_1.Fill( evt.pt() )
        metChargeHist_1.Fill( evt.charge() )
    
    for evt in tauH.product():
        print "Tau       pT %8.2f " % evt.pt()
        iTauCount_1+= 1
        tauPtHist_1.Fill( evt.pt() )
        tauChargeHist_1.Fill( evt.charge() )

    for evt in jetH.product():
        print "Jet  pT %8.2f " % evt.pt()
        jetPtHist_1.Fill( evt.pt() )
        jetChargeHist_1.Fill( evt.charge() )

#    for evt in eleH.product():
#        print "Electron  pT %8.2f " % evt.pt()
#        elecPtHist_1.Fill( evt.pt() )
#        elecChargeHist_1.Fill( evt.charge() )
    print "-----"

muonH, muonN = Handle("std::vector<pat::Muon>"), "MEMProducer"
metH, metN = Handle("std::vector<pat::MET>"), "MEMProducer"
tauH, tauN = Handle("std::vector<pat::Tau>"), "MEMProducer"
jetH, jetN = Handle("std::vector<pat::Jet>"), "MEMProducer"
#eleH, eleN = Handle("std::vector<pat::Electron>"), "MEMProducer"
chiH, chiN = Handle("std::vector<reco::MEMResult>"), "MEMProducer"

events = Events([],options=options)
iTauCount_2 = 0
iMuonCount_2 = 0
print "\nboucle enumerate 2"
for i,event in enumerate(events):
    print "Event : ", i
    event.getByLabel(muonN, muonH)
    event.getByLabel(metN, metH)
    event.getByLabel(tauN, tauH)
    event.getByLabel(jetN, jetH)
#    event.getByLabel(eleN, eleH)
    event.getByLabel(chiN, chiH)
    
    for evt in muonH.product():
        print "Muon      pT %8.2f " % evt.pt()
        iMuonCount_2 += 1
        muonPtHist_2.Fill( evt.pt() )
        muonChargeHist_2.Fill( evt.charge() )

    for evt in metH.product():
        print "MET       pT %8.2f " % evt.pt()
        metPtHist_2.Fill( evt.pt() )
        metChargeHist_2.Fill( evt.charge() )

    for evt in tauH.product():
        print "Tau       pT %8.2f " % evt.pt()
        iTauCount_2 += 1
        tauPtHist_2.Fill( evt.pt() )
        tauChargeHist_2.Fill( evt.charge() )

    for evt in jetH.product():
        print "Jet  pT %8.2f " % evt.pt()
        jetPtHist_2.Fill( evt.pt() )
        jetChargeHist_2.Fill( evt.charge() )

#    for evt in eleH.product():
#        print "Electron  pT %8.2f " % evt.pt()
#        elecPtHist_2.Fill( evt.pt() )
#        elecChargeHist_2.Fill( evt.charge() )

    for evt in chiH.product():
        print "Chi nb of daughters : %i " % evt.numberOfDaughters()
        print "Chi muon  pT %8.2f " % evt.daughter(0).pt()
        print "Chi tau   pT %8.2f " % evt.daughter(1).pt()
        print "Chi jet 1 pT %8.2f " % evt.daughter(2).pt()
        print "Chi jet 2 pT %8.2f " % evt.daughter(3).pt()
        print "Chi met   pT %8.2f " % evt.daughter(4).pt()

    print "-----"

print "Tau - count 1 : ", iTauCount_1, " - count 2 : ", iTauCount_2
print "Muon - count 1 : ", iMuonCount_1, " - count 2 : ", iMuonCount_2
  
f_rel = ROOT.TFile("2016.11.24_17.12.53/myOutputFile.ID_0.root")
f_rel.GetListOfKeys().Print()
f_rel.ls()
hname="Events"
tt=f_rel.FindKey(hname)
if (tt):
    print "OK : ", tt.ls()
    print "OK : ", tt.Print()

createPicture(muonPtHist_1, muonPtHist_2, "muonPtHist_py.png")
createPicture(muonChargeHist_1, muonChargeHist_2, "muonChargeHist_py.png")
createPicture(metPtHist_1, metPtHist_2, "metPtHist_py.png")
createPicture(metChargeHist_1, metChargeHist_2, "metChargeHist_py.png")
createPicture_2(tauPtHist_1, tauPtHist_2, "tauPtHist_py.png")
createPicture_2(tauChargeHist_1, tauChargeHist_2, "tauChargeHist_py.png")
createPicture_2(jetPtHist_1, jetPtHist_2, "jetPtHist_py.png")
createPicture_2(jetChargeHist_1, jetChargeHist_2, "jetChargeHist_py.png")
#createPicture(elecPtHist_1, elecPtHist_2, "elecPtHist_py.png")
#createPicture(elecChargeHist_1, elecChargeHist_2, "elecChargeHist_py.png")

print "Fin."
