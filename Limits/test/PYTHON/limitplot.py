#!/usr/bin/env python
import sys, os
import ROOT
from CombineHarvester.CombineTools.plotting import *
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(ROOT.kTRUE)

# Style and pads
ModTDRStyle()
graphs = StandardLimitsFromJSONFile(sys.argv[1])
name = os.path.basename(sys.argv[1])
sname = name.replace(".json","")
print sname
canv = ROOT.TCanvas(sname, sname)
pads = OnePad()

# Get limit TGraphs as a dictionary

# Create an empty TH1 from the first TGraph to serve as the pad axis and frame
axis = CreateAxisHist(graphs.values()[0])
axis.GetXaxis().SetTitle('m_{res} (GeV)')
axis.GetYaxis().SetTitle('95% CL limit')
pads[0].cd()
pads[0].SetLogy()
axis.Draw('axis')

# Create a legend in the top left
#legend = PositionedLegend(0.3, 0.2, 3, 0.015)
legend = PositionedLegend(0.3, 0.2, 3, 0.015)

# Set the standard green and yellow colors and draw
StyleLimitBand(graphs)
DrawLimitBand(pads[0], graphs, legend=legend)
legend.Draw()

# Re-draw the frame and tick marks
pads[0].RedrawAxis()
pads[0].GetFrame().Draw()

# Adjust the y-axis range such that the maximum graph value sits 25% below
# the top of the frame. Fix the minimum to zero.
FixBothRanges(pads[0], 0.1, 0, GetPadYMax(pads[0]), 0.25)

if "0p01" in sys.argv[1]:
    text = "width = 0.01%"
if "5" in sys.argv[1]:
    text = "width = 5%"
text2=None
if "mu" in  sname or "el" in sname:
    ch = sname.split("_")[-1][:2]
    year = int(sname.split("_")[-1][2:])
    ch2 = "Electron" if ch=="el" else "Muon"
    text2 = "%s Channel %i" %(ch2,year)
    print text2
l = ROOT.TLatex()
l.SetNDC()
l.SetTextSize(0.04)
l.SetTextFont(72)
#l.DrawLatex(0.6,.65, "W#gamma#rightarrow l#gamma#nu")
#l.DrawLatex(0.6,.6, text)
if text2: l.DrawLatex(0.6,.25, text2)

hmax = axis.GetMaximum()
#hmin = axis.GetMinimum()
hmin = 0.1
def oneline(xval):
    ln = ROOT.TLine(xval,hmax,xval,hmin)
    ln.SetLineStyle(3)
    ln.SetLineWidth(2)
    ln.Draw()
    return ln
ln = oneline(625)
#ln2 = oneline(425)

# Standard CMS logo
DrawCMSLogo(pads[0], 'CMS', 'Simulation Preliminary', 11, 0.200, 0.035, 1.2, '', 0.8)



canv.Print('.pdf')
canv.Print('.png')
