import FWCore.ParameterSet.Config as cms
process = cms.Process("Demo")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) ) 
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load("HeterogeneousCore.CUDAServices.CUDAService_cfi")
process.GlobalTag.globaltag = 'GR_P_V42_AN3::All'
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T21', '')
process.MessageLogger.cerr.FwkReport.reportEvery = 1
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(10)) 
process.source = cms.Source("PoolSource",
# replace 'myfile.root' with the source file you want to use 
fileNames = cms.untracked.vstring(
#'file:/afs/cern.ch/user/g/gpizzati/CMSSW_12_2_0_pre3/src/step2.root' 
'/store/relval/CMSSW_12_3_0_pre5/RelValTTbar_14TeV/GEN-SIM-RECO/123X_mcRun4_realistic_v4_2026D88noPU-v1/10000/7cbdf153-c397-410e-8af5-dc3905565ced.root'
)
)

process.demo = cms.EDProducer(
    'gpuTest',
    TrackLabel = cms.InputTag("generalTracks"),
    beamSpotLabel = cms.InputTag("offlineBeamSpot"),

)
process.out = cms.OutputModule("PoolOutputModule", 
    fileName= cms.untracked.string("test3.root")
)
"""
process.TFileService = cms.Service("TFileService",
fileName = cms.string('outfile.root')
)
"""
process.p = cms.Path(process.demo)
process.ep = cms.EndPath(process.out)
