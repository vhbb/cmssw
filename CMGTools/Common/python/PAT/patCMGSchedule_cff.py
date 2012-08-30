import FWCore.ParameterSet.Config as cms
from CMGTools.Common.Tools.cmsswRelease import *

#this one
def getSchedule(process, runOnMC, runOnFastSim):

    result = cms.Schedule(
        process.p,
        process.EcalDeadCellTriggerPrimitiveFilterPath,
        process.hcalLaserEventFilterPath,
        process.trackingFailureFilterPath,
        process.primaryVertexFilterPath,
        process.noscrapingFilterPath,
        # doesn't work in 44X 
        # process.eeBadScFilterPath,
        process.metNoiseCleaningPath
        )
    if runOnMC:
        result.append(process.totalKinematicsFilterPath)
    if not( runOnFastSim ):
        result.append(process.CSCTightHaloFilterPath)
        result.append(process.HBHENoiseFilterPath)
    else:
        process.metNoiseCleaningPath.remove(process.CSCTightHaloFilter)
        process.metNoiseCleaningPath.remove(process.HBHENoiseFilter)
    return result

