
# coding: utf-8

# In[1]:

import ROOT as r


# In[2]:

#r.gSystem.Load("RooDCB.C")
r.gROOT.ProcessLine(".L RooDCB.C")


# In[3]:

import cPickle
import numpy as np
import pandas as pd
import time


# In[4]:

FIT_LOW, FIT_HIGH = 105e3, 160e3
#FIT_LOW, FIT_HIGH = 120e3, 130e3

CAT = 2


# In[5]:

def draw():
    r.gROOT.FindObject('c1').Draw()
    r.gROOT.FindObject('c1').Update()

def logy(dolog=True):
    r.gROOT.FindObject('c1').SetLogy(dolog)


# In[6]:

myy = r.RooRealVar("myy","myy",105e3,160e3)


# # Load in m_yy histograms from Xhh_m350_yybb (h014)

# In[7]:

f = r.TFile("histos_SM.root")
keys = [k.GetName() for k in f.GetListOfKeys()]
histos_all = [{},{},{}]
datasets_all = [{},{},{}]
for k0 in keys:
    if not k0.startswith('h_'): continue
    k = k0[len('h_'):]
    cat = int(k.split('_')[-1])
    v = '_'.join(k.split('_')[:-1])
    h = f.Get(k0)
    histos_all[cat][v] = h
    
    ds = r.RooDataHist('d_'+k,'d_'+k,r.RooArgList(myy),h)
    datasets_all[cat][v] = ds


# In[8]:

m0 = r.RooRealVar("m0","m0",120e3,130e3)
sigma0 = r.RooRealVar("sigma0","sigma0",0,5e3)
#alphaHi = r.RooRealVar("alphaHi","alphaHi",1.886)#0,20)
alphaHi = r.RooRealVar("alphaHi","alphaHi",0,20)
nHi = r.RooRealVar("nHi","nHi",2.42)#0,20)
#nHi = r.RooRealVar("nHi","nHi",0,20)
#alphaLo = r.RooRealVar("alphaLo","alphaLo",1.591)#0,20)
alphaLo = r.RooRealVar("alphaLo","alphaLo",0,20)
nLo = r.RooRealVar("nLo","nLo",4.23)#0,20)
#nLo = r.RooRealVar("nLo","nLo",0,20)

dcb = r.RooDCB("dcb","dcb",myy,m0,sigma0,alphaHi,nHi,alphaLo,nLo)


# In[9]:

r.RooMsgService.instance().setGlobalKillBelow(r.RooFit.FATAL)


# In[10]:

ds = datasets_all[2]['nominal']
fit_result = dcb.fitTo(ds,r.RooFit.Verbose(False))


# In[11]:

frame = myy.frame()
ds.plotOn(frame)
dcb.plotOn(frame)
frame.SetMinimum(1e-3)
frame.Draw()
logy(True)
draw()


# In[12]:

histos = histos_all[CAT]
all_variations = sorted(histos.keys())
up_variations = [v for v in all_variations if v.lower().endswith('up') and not 'JER' in v]
down_variations = [v for v in all_variations if v.lower().endswith('down')]
other_variations = list(set(all_variations).difference(up_variations+down_variations))

print("{} variations: {} up / {} down / {} other".format(len(all_variations), len(up_variations), len(down_variations), len(other_variations)))


# In[ ]:

ds_nom = datasets_all[CAT]['nominal']
dcb.fitTo(ds_nom)
m0_nom = m0.getVal()
sigma0_nom = sigma0.getVal()
integral_nom = histos_all[CAT]['nominal'].Integral()


# In[ ]:

fit_params = []
fits_up = []
hists_up = []
fits_down = []
hists_down = []
up_bases = []
down_bases = []
for v in sorted(datasets_all[CAT].keys()):
    if v == 'nominal': continue
    
    ds = datasets_all[CAT][v]
    #h.Fit(f_dg, 'Q', '', FIT_LOW, FIT_HIGH)
    fit_result = dcb.fitTo(ds, r.RooFit.Extended(True), r.RooFit.Save())
    status = fit_result.status()

    print "Got status = ", status
    #if status != 0: break
    if status != 0:
        frame = myy.frame()
        ds.plotOn(frame)
        dcb.plotOn(frame)
        frame.SetMinimum(1e-4)
        frame.Draw()
        logy(True)
        draw()
        fit_result.Print()
        raw_input("Press enter to continue.")
    
    integral = histos_all[CAT][v].Integral()
    
    sign = 0
    if v in up_variations:
        sign = 1
    if v in down_variations:
        sign = -1
        
    #result = (m0,c1,sigma1,c2,sigma2,sign,v)
    #result = (m0.getVal(),c1,sigma1,c2,sigma2,integral,sign,v)
    result = (m0.getVal(),sigma0.getVal(),integral,sign,v,status)
    fit_params.append(result)
    
    if v == 'nominal':
        nominal_fit = np.array(result)
        continue
    
    if v in up_variations:
        fits_up.append(result)
        hists_up.append(h)
        up_bases.append(v.split('__')[0])
    if v in down_variations:
        fits_down.append(result)
        hists_down.append(h)
        down_bases.append(v.split('__')[0])

cols = ['m0','sigma0','integral','sign','variation','status']
fit_params = pd.DataFrame(fit_params, columns=cols)


# In[ ]:

# calculate fractional errors w.r.t nominal
fit_params = fit_params.assign(dm0=lambda x: 100.*(x.m0-m0_nom)/m0_nom,
                               dsigma1=lambda x: 100.*(x.sigma0-sigma0_nom)/sigma0_nom,
                               dnorm=lambda x: 100.*(x.integral-integral_nom)/integral_nom)


# In[ ]:

# flag variations w/ more than 0.1% effect
fit_params = fit_params.assign(m0_flag=lambda x: np.abs(x.dm0)>0.1,
                               sigma1_flag=lambda x: np.abs(x.dsigma1)>0.1,
                               norm_flag=lambda x: np.abs(x.dnorm)>0.1)


# # dominant systematics on signal width (%)

print "dominant systematics on signal width (%)"
print "="*40
print fit_params[fit_params.sigma1_flag][['variation','dsigma1','status']]
print "="*40
print


# # dominant systematics on signal mean (%)

print "dominant systematics on signal mean (%)"
print "="*40
print fit_params[fit_params.m0_flag][['variation','dm0','status']]
print "="*40
print


print "dominant systematics on signal yield (%)"
print "="*40
print fit_params[fit_params.norm_flag][['variation','dnorm','status']]
print "="*40
print

#print fit_params[['variation','dnorm']]
#print fit_params[['variation','status']]

print "Nonzeros statuses: ", np.sum(fit_params.status != 0)
