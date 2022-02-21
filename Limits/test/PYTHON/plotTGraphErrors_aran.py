#!/usr/bin/env python
# Aran Garcia-Bellido (Dec 7, 2010)
#
# This script will read in one or more text files with columns of x y numbers and make a scatter plot with them.
# If several files are given in the input, the xy lines are overlaid. 
# There is flexibility in the format of the input files: there can be more than one space between columns, 
# and there can be several columns. 2 columns is just x y; 3 columns means x y error_y; 
# and 4 columns means x y error_x error_y; 
# Additionally, the files can have comments in lines starting with #
#
import sys
import string
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from ROOT import TCanvas, gInterpreter, gROOT, TLegend, TGraph, TGraphErrors, TMultiGraph, TH1F, nullptr

gROOT.SetBatch(True)

def main():
	parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
	parser.add_argument("-n", "--name", dest="name", type=str, required=True, help="name of pdf to be saved")
	parser.add_argument("-f", "--files", dest="files", type=str, required=True, nargs='+', help="list of all the files to plot, each w/ columns x y [err]")
	parser.add_argument("-l", "--doLimit", dest="doLimit", default=False, action="store_true", help="limit mode (use combined regions)")
	parser.add_argument("-x", "--extra", dest="extra", type=str, default="", help="extra text for legend")
	args = parser.parse_args()

	# Specify stuff to make pretty plots (color, marker and name for the legend)
	color_list=[4,2,1,6,7]
	marker_list=[20,21,22,23,24]
	
	if args.doLimit:
		name_list=["cut", "bdt"]
	else:
		name_list=["lowCut", "lowSVJ2","highCut","highSVJ2"]
	files_dim = len(args.files)
	print " Will process %i files " % files_dim 
	
	# loop over files:
	allGraphs = TMultiGraph()
	leg = TLegend(0.55, 0.70, 0.89, 0.89) # upper right
	leg.SetBorderSize(0)
	leg.SetTextSize(0.04)
	leg.SetFillColor(0)
	i=0
	if len(args.extra)>0:
		leg.AddEntry(nullptr,args.extra,"")
	for f in args.files:
		g=Get_Graph_from_file(f)
		g.SetLineColor(color_list[i])
		g.SetMarkerColor(color_list[i])
		g.SetMarkerStyle(marker_list[i])
		g.SetFillColor(0)
		g.SetLineWidth(2)
		allGraphs.Add(g,"lp")
		leg.AddEntry(g,name_list[i],"lp")
		i=i+1
	if "mean" in args.name:
		allGraphs.SetMaximum(1.5)
		allGraphs.SetMinimum(-1.5)
	elif "stdev" in args.name:
		allGraphs.SetMaximum(2.5)
		allGraphs.SetMinimum(-0.5)
	elif "chi2" in args.name:
		allGraphs.SetMaximum(5.0)
		allGraphs.SetMinimum(0.0)
	else:
		pass
		
	gROOT.SetStyle('Plain')
	c1 = TCanvas("c1", "ReadSeveralFilesScatterPlotErrors")
	c1.SetGrid(1,1)
	c1.GetFrame().SetFillColor(0)
   	c1.GetFrame().SetBorderSize(0)
   	#c1.SetLogy()
   
	#allGraphs.SetMaximum(ymax)
	#allGraphs.SetMinimum(0.0)
	allGraphs.Draw("0AL")
	leg.Draw()
	if "Sig1" in args.name:
		ring = 1
	elif "Sig0" in args.name:
		ring = 0
	elif "SigM" in args.name:
		ring = "median r_{exp}"
	elif "Sig3M" in args.name:
		ring = "3 * median r_{exp}"
	if ("GenMainFitMain" in args.name) or ("GenAltFitAlt" in args.name):
		extraTitle = ", Closure"
	else:
		extraTitle = ", Bias"
	if "GenMain" in args.name:
		extraTitle += ", GenMain"
	else:
		extraTitle += ", GenAlt"
	if "mean" in args.name:
		allGraphs.SetTitle("Mean of Gaus Fit to Pull Distribution, r = {}".format(ring) + extraTitle)
		allGraphs.GetYaxis().SetTitle("Mean")
	elif "stdev" in args.name:
		allGraphs.SetTitle("StdDev of Gaus Fit to Pull Distribution, r = {}".format(ring) + extraTitle)
		allGraphs.GetYaxis().SetTitle("StdDev")
	elif "chi2" in args.name:
		allGraphs.SetTitle("#chi^2 of Gaus Fit to Pull Distribution, r = {}".format(ring) + extraTitle)
		allGraphs.GetYaxis().SetTitle("#chi^2")
		
	allGraphs.GetXaxis().SetAxisColor(17)
	allGraphs.GetYaxis().SetAxisColor(17)
	if "varyZ" in args.name:
		allGraphs.GetXaxis().SetTitle("m_{Z'} [GeV]")
	if "varyD" in args.name:
		allGraphs.GetXaxis().SetTitle("m_{D} [GeV]")
	if "varyR" in args.name:
		allGraphs.GetXaxis().SetTitle("r_{inv}")
	if "varyA" in args.name:
		allGraphs.GetXaxis().SetTitle("#alpha_{D}")
	#allGraphs.GetYaxis().SetTitle("#sigma(p#bar{p}#rightarrow#phi) #times Br(#phi#rightarrow#tau#bar{#tau})) [pb]")
	c1.RedrawAxis()
	c1.Update()
	c1.SaveAs(args.name+".pdf")
	
def Get_Graph_from_file(filename="table.dat"):
	from array import array
	
	xx	 = array('f')
	yy	 = array('f')
	err_xx = array('f')
	err_yy = array('f')

	cdim = 0 # prevents error later on
	
	for line in open(filename):
		line.replace('\n',' ').replace('\r',' ') # remove carriage returns
		if '#' in line: continue # skip comment lines on the file
		columns = line.split()   # each line is read into a list called columns
		cdim=len(columns)		# how many columns does this file have?
		if cdim < 2: continue	# this prevents crashing if there is an empty line at the end
		xx.append(eval(columns[0]))
		yy.append(eval(columns[1]))
		if cdim == 3: # x y erry
			err_yy.append(float(columns[2]))
			err_xx.append(0.0)
		elif cdim == 4: # x y errx erry
			err_xx.append(eval(columns[2]))	
			err_yy.append(float(columns[3]))
		else:
			err_xx.append(0.0)	
			err_yy.append(0.0)
	
	#print xx,yy,err_xx,err_yy
	nlines=len(xx)
	if cdim == 0:
		print "{} was empty or DNE".format(filename)
		return TGraphErrors()
	print " File %-20s read in with %i rows and %i columns " %(filename,nlines,cdim)	
	
	gr = TGraphErrors(nlines,xx,yy,err_xx,err_yy)
	gr.Sort()
	#gr.SetTitle("TGraphErrors Example")
	#gr.GetHistogram().SetXTitle("m_{H} [GeV]")
	#gr.GetHistogram().SetYTitle("#sigma(#p#bar{p}#rightarrow#phi) #times Br(#phi#rightarrow#tau#bar{tau})) [pb]")
	return gr
		
main()
