import sys

fname = sys.argv[1]
wsname = sys.argv[2]
pdfname1 = sys.argv[3]
pdfname2 = sys.argv[4]
highCut3 = len(sys.argv)>5 and "highCut" in fname

sys.argv.append('-b')
import ROOT as r
r.gROOT.ProcessLine("""class RooGenericPdf2 : public RooGenericPdf {
public:
RooListProxy actualVars() { return _actualVars; }
TString formExpr() { return _formExpr; }
static RooListProxy getActualVars(RooGenericPdf* pdf){ return ((RooGenericPdf2*)pdf)->actualVars(); }
static TString getFormExpr(RooGenericPdf* pdf){ return ((RooGenericPdf2*)pdf)->formExpr(); }
};""")
xmin = 1500
xmax = 8000
nbins = 65
file = r.TFile.Open(fname)
ws = file.Get(wsname)
pdf_main = ws.pdf(pdfname1)
args_main = r.RooGenericPdf2.getActualVars(pdf_main)
expr_main = str(r.RooGenericPdf2.getFormExpr(pdf_main))
pdf_alt = ws.pdf(pdfname2)
args_alt = r.RooGenericPdf2.getActualVars(pdf_alt)
expr_alt = str(r.RooGenericPdf2.getFormExpr(pdf_alt))

def getParValues(ws,args):
    iter = args.createIterator()
    var = iter.Next()
    ctr = 0
    vals = []
    while var:
        if ctr>0: vals.append(ws.var(var.GetName()).getVal())
        ctr += 1
        var = iter.Next()
    return vals

def fixExpr(expr):
    expr = expr.replace("@0","x")
    expr2 = expr[:]
    expr2_ind = 0
    for ind in range(len(expr)):
        if expr[ind]=='@':
            expr2 = expr2[:expr2_ind] + "[{}]".format(expr[ind+1]) + expr2[expr2_ind+2:]
            expr2_ind += 1
        expr2_ind += 1
    expr2 = "[0]*({})".format(expr2)
    return expr2

expr_main = fixExpr(expr_main)
f_main = r.TF1("f_main",expr_main,xmin,xmax)
v_main = getParValues(ws,args_main)
print "v_main = {}".format(v_main)
f_main.SetParameter(0,1)
for iv,v in enumerate(v_main):
    f_main.SetParameter(iv+1,v)
data_yield = ws.data("data_obs").sumEntries()*100
f_main.SetParameter(0, data_yield/f_main.Integral(xmin,xmax))
hist = r.TH1F("Bkg","",nbins,xmin,xmax)
hist.Add(f_main)
hist.GetXaxis().SetTitle("M_{T} [GeV]")
hist.GetYaxis().SetTitle("number of events")
for b in range(nbins): hist.SetBinError(b+1, r.TMath.Sqrt(hist.GetBinContent(b+1)))

mopt = r.Math.MinimizerOptions()
mopt.SetMaxFunctionCalls(100000)
fopt = r.Foption_t()
r.Fit.FitOptionsMake(r.Fit.kHistogram,"QN0",fopt)
ropt = r.Fit.DataRange(xmin,xmax)

expr_alt = fixExpr(expr_alt)
f_orig = r.TF1("f_orig",expr_alt,xmin,xmax)
v_alt = getParValues(ws,args_alt)
print "v_alt = {}".format(v_alt)
f_orig.SetParameter(0,1)
for iv,v in enumerate(v_alt):
    f_orig.SetParameter(iv+1,v)
data_yield = ws.data("data_obs").sumEntries()*100
f_orig.SetParameter(0, data_yield/f_orig.Integral(xmin,xmax))
f_orig.SetLineColor(r.kBlue)

if highCut3:
    expr_alt = "[0]*(exp([1]*(x/13000)) * pow(x/13000,[2]*(1+[3]*log(x/13000))))"
    v_alt.append(1)

def varyAll(pos,paramlist,val,tups):
    vals = paramlist[pos]
    for v in vals:
        tmp = val[:]+[v]
        # check if last param
        if pos+1==len(paramlist):
            tups.add(tuple(tmp))
        else:
            varyAll(pos+1,paramlist,tmp,tups)

# norm only gets positive values
init_vals = [
    [0.1,1,10],
]
for i in range(len(v_alt)):
    init_vals.append([-10,-1,-0.1,0.1,1,10])
# set to accumulate all init combos
initset = set()
varyAll(0,init_vals,[],initset)
inits = list(sorted(initset))

def setInit(fn,init):
    for p,par in enumerate(init):
        fn.SetParameter(p,par)
        if p>0: fn.SetParLimits(p,-100,100)    

best_init = None
best_chi2 = 1e100
for i,init in enumerate(inits):
    f_alt = r.TF1("f_alt_{}".format(i),expr_alt,xmin,xmax)
    setInit(f_alt,init)
    result = r.Fit.FitObject(hist,f_alt,fopt,mopt,"",ropt)
    chi2 = f_alt.GetChisquare()
    if chi2 < best_chi2:
        best_chi2 = chi2
        best_init = init

# now refit and keep the best one
print best_init
f_alt = r.TF1("f_alt",expr_alt,xmin,xmax)
setInit(f_alt,best_init)
fopt2 = r.Foption_t()
r.Fit.FitOptionsMake(r.Fit.kHistogram,"0",fopt2)
status = r.Fit.FitObject(hist,f_alt,fopt2,mopt,"",ropt)

can1 = r.TCanvas("c1")
hist.Draw()
f_orig.Draw("same")
f_alt.Draw("same")
can1.SetLogy()
can1.Update()
can1.Print("fit_best__{}__{}{}.png".format(pdfname1,pdfname2,"_more" if highCut3 else ""),"png")

def makeRes(hist, hname, fn, color):
    hres = hist.Clone(hname)
    hres.GetFunction("f_alt").Delete()
    hres.Add(fn,-1)
    ymax = max(abs(hres.GetMaximum()),abs(hres.GetMinimum()))
    hres.SetLineColor(color)
    hres.SetMarkerColor(color)
    hres.GetXaxis().SetTitle("M_{T} [GeV]")
    hres.GetYaxis().SetTitle("residual (data-fit)")
    return hres,ymax
    
hres_orig, ymax_orig = makeRes(hist,"hres_orig",f_orig,r.kBlue)
hres_alt, ymax_alt = makeRes(hist,"hres_alt",f_alt,r.kRed)
ymax = max(ymax_orig, ymax_alt)

can2 = r.TCanvas("c2")
hres_orig.GetYaxis().SetRangeUser(-ymax,ymax)
hres_orig.Draw()
hres_alt.Draw("same")
can2.Print("res_best__{}__{}{}.png".format(pdfname1,pdfname2,"_more" if highCut3 else ""),"png")

