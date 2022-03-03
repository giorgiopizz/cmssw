import FWCore.ParameterSet.Config as cms
process = cms.Process("Demo")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) ) 
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = 'GR_P_V42_AN3::All'
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T21', '')
process.MessageLogger.cerr.FwkReport.reportEvery = 10
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1)) 
process.source = cms.Source("PoolSource",
# replace 'myfile.root' with the source file you want to use 
fileNames = cms.untracked.vstring(
'file:/afs/cern.ch/user/g/gpizzati/CMSSW_12_2_0_pre3/src/step3.root' 
)
)

process.demo = cms.EDAnalyzer('MyTransientTrackBuilderTest')
"""
process.TFileService = cms.Service("TFileService",
fileName = cms.string('outfile.root')
)
"""
process.p = cms.Path(process.demo)

