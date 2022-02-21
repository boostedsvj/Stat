import ROOT as rt
#rt.gROOT.SetBatch(1)
#rt.gROOT.SetOptStats(0)
rt.gStyle.SetOptFit(1111)

# this should be run in the Stats/Limits/test folder.
rt.gSystem.Load("libHiggsAnalysisCombinedLimit.so")
cardsDirect = "cards_06Jan/"
signalDirect = "SVJ_mZprime3000_mDark20_rinv03_alphapeak/"
for toyStr in [""]:#, "DijetStr0"]:#, "modExpStr0", "modExpStr1", "ExpStr0", "ExpStr1"]:
	for region in ["lowSVJ1","lowSVJ2","highSVJ1","highSVJ2"]:
		_file = rt.TFile(cardsDirect+signalDirect+region+"/fitDiagnostics"+toyStr+".root","read")
		tree = _file.Get("tree_fit_sb")
		tree.Draw("(r-1)/rErr>>h(20,-1,1)")
		tree.Draw("(r-1)/rErr>>hfs(20,-1,1)","fit_status>=0")
		h = rt.gDirectory.Get("h")
		hfs = rt.gDirectory.Get("hfs")
		print("Fits in full:{}, fits not failed: {}".format(h.GetEntries(), hfs.GetEntries()))
		h.SetDirectory(0)
		h.SetName("hist_"+region)
		h.Fit("gaus","Q")
		hfs.SetDirectory(0)
		hfs.SetName("hist_"+region+"_fs0")
		hfs.Fit("gaus","Q")
		c = rt.TCanvas("c","c",900,600)
		h.Draw()
		c.SaveAs(cardsDirect+signalDirect+"plots/"+h.GetName()+"_"+toyStr+".png")
		hfs.Draw()
		c.SaveAs(cardsDirect+signalDirect+"plots/"+hfs.GetName()+"_"+toyStr+".png")
		_file.Close()
		
		
