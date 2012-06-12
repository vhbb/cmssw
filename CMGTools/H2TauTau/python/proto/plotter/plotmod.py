import copy
from CMGTools.H2TauTau.proto.plotter.H2TauTauDataMC import H2TauTauDataMC
from CMGTools.RootTools.Style import *
from ROOT import kPink

def buildPlot( var, anaDir,
               comps, weights, nbins, xmin, xmax,
               cut, weight,
               embed, treeName ):
    pl = H2TauTauDataMC(var, anaDir,
                        comps, weights, nbins, xmin, xmax,
                        str(cut), weight,
                        embed, treeName )
    return pl


def hist( var, anaDir,
          comp, weights, nbins, xmin, xmax,
          cut, weight,
          embed, treeName=None ):
    pl = buildPlot( var, anaDir,
                    {comp.name:comp}, weights, nbins, xmin, xmax,
                    cut, weight,
                    embed, treeName )
    histo = copy.deepcopy( pl.Hist(comp.name) )    
    return histo


def shape( var, anaDir,
           comp, weights, nbins, xmin, xmax,
           cut, weight,
           embed, treeName=None):
    shape = hist( var, anaDir,
                  comp, weights, nbins, xmin, xmax,
                  cut, weight,
                  embed, treeName )
    shape.Normalize()
    return shape


def shape_and_yield( var, anaDir,
                     comp, weights, nbins, xmin, xmax,
                     cut, weight,
                     embed, treeName=None ):
    shape = hist( var, anaDir,
                  comp, weights, nbins, xmin, xmax,
                  cut, weight,
                  embed, treeName )
    yi = shape.Integral()
    shape.Normalize()
    return shape, yi

    
def addQCD( plot, dataName ):
    # import pdb; pdb.set_trace()
    plotWithQCD = copy.deepcopy( plot )
    # QCD_data = data - DY - TTbar - W
    qcd = copy.deepcopy(plotWithQCD.Hist(dataName))
    qcd.Add(plotWithQCD.Hist('Ztt'), -1)
    try:
        dyJetsFakes = plot.Hist('Ztt_Fakes')
        qcd.Add(dyJetsFakes, -1)
    except:
        print 'cannot find Ztt_Fakes'
        print plot
        pass    
    qcd.Add(plotWithQCD.Hist('TTJets'), -1)
    qcd.Add(plotWithQCD.Hist('WJets'), -1)
    # adding the QCD data-driven estimation to the  plot
    plotWithQCD.AddHistogram( 'QCD', qcd.weighted, 888)
    plotWithQCD.Hist('QCD').stack = True
    plotWithQCD.Hist('QCD').SetStyle( sHTT_QCD )
    return plotWithQCD


def getQCD( plotSS, plotOS, dataName ):

    # use SS data as a control region
    # to get the expected QCD shape and yield
    plotSSWithQCD = addQCD( plotSS, dataName )

    # extrapolate the expected QCD shape and yield to the
    # signal region

    plotOSWithQCD = copy.deepcopy( plotOS )

    qcdOS = copy.deepcopy( plotSSWithQCD.Hist('QCD') )
    qcdOS.RemoveNegativeValues()
    qcdOS.Scale( 1.11 )

    plotOSWithQCD.AddHistogram('QCD', qcdOS.weighted, 1030)
    plotOSWithQCD.Hist('QCD').layer=1.5

    return plotSSWithQCD, plotOSWithQCD


def groupEWK( plot ):
    wjets = plot.Hist('WJets')
    ewk = copy.deepcopy( wjets )
    dyfakes = plot.Hist('Ztt_Fakes')
    dyfakes.on = False
    wjets.on = False
    ewk.Add( dyfakes )
    plot.AddHistogram('EWK', ewk.weighted, dyfakes.layer, 'EWK') 



def fW_inclusive(mtplot, dataName, xmin, xmax):
    
    # WJets_data = data - DY - TTbar
    wjet = copy.deepcopy(mtplot.Hist(dataName))
    wjet.Add(mtplot.Hist('Ztt'), -1)
    try:
        dyJetsFakes = mtplot.Hist('Ztt_Fakes')
        wjet.Add(mtplot.Hist('Ztt_Fakes'), -1)
    except:
        pass
    # FIXME
    wjet.Add(mtplot.Hist('TTJets'), -1)

    # adding the WJets_data estimation to the stack
    mtplot.AddHistogram( 'Data - DY - TT', wjet.weighted, 1010)
    mtplot.Hist('Data - DY - TT').stack = False
    # with a nice pink color
    pink = kPink+7
    sPinkHollow = Style( lineColor=pink, markerColor=pink, markerStyle=4)
    mtplot.Hist('Data - DY - TT').SetStyle( sPinkHollow )

    # determine scaling factor for the WJet MC
    mtmin, mtmax = xmin, xmax
    # scale = WJets_data / WJets 
    scale_WJets = mtplot.Hist('Data - DY - TT').Integral(True, mtmin, mtmax) \
                  / mtplot.Hist('WJets').Integral(True, mtmin, mtmax)
    # apply this additional scaling factor to the WJet component 
    mtplot.Hist('WJets').Scale(scale_WJets)

    # hide the WJets_data component from the mtplot. can be set to True interactively
    mtplot.Hist('Data - DY - TT').on = True

    mtplot.Hist('WJets').layer = -999999
    return scale_WJets



def plot_W_inclusive(var, anaDir,
                     comps, weights, nbins, xmin, xmax,
                     cut, weight,
                     embed):

    # get WJet scaling factor for same sign
    print 'extracting SS WJets inclusive data/MC factor'
    var = 'mt'
    sscut = 'isSignal && mt>{mtcut} && diTau_charge!=0'.format(mtcut=xmin)
    oscut = 'isSignal && mt>{mtcut} && diTau_charge==0'.format(mtcut=xmin)
    mtSS = H2TauTauDataMC(var, anaDir, comps, weights,
                          nbins, xmin, xmax,
                          cut = sscut, weight=weight,
                          embed=embed)
    # replaceWJetShape( mtSS, var, sscut)
    # import pdb; pdb.set_trace()
    fW_inclusive_SS = fW_inclusive( mtSS, 'Data', xmin, xmax)
    # get WJet scaling factor for opposite sign
    print 'extracting OS WJets inclusive data/MC factor'
    mtOS = H2TauTauDataMC(var, anaDir, comps, weights,
                          nbins, xmin, xmax, 
                          cut = oscut, weight=weight,
                          embed=embed)
    # replaceWJetShape( mtOS, var, oscut)
    fW_inclusive_OS = fW_inclusive( mtOS, 'Data', xmin, xmax)
    print 'fW_inclusive_SS=',fW_inclusive_SS,'fW_inclusive_OS=',fW_inclusive_OS
    return fW_inclusive_SS, fW_inclusive_OS, mtSS, mtOS

