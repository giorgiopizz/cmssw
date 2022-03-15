import FWCore.ParameterSet.Config as cms
from RecoVertex.PrimaryVertexProducer.OfflinePrimaryVerticesCUDA_cfi import offlinePrimaryVertices
process = cms.Process("Demo")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) ) 
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
#process.load('HeterogeneousCore.CUDACore.ProcessAcceleratorCUDA_cfi')
process.load("HeterogeneousCore.CUDAServices.CUDAService_cfi")

#process.load( "HLTrigger.Timer.FastTimerService_cfi" )

process.GlobalTag.globaltag = 'GR_P_V42_AN3::All'
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T21', '')
process.MessageLogger.cerr.FwkReport.reportEvery = 1
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(30)) 
"""
process.IgProfService = cms.Service("IgProfService",
  reportFirstEvent            = cms.untracked.int32(0),
  reportEventInterval         = cms.untracked.int32(25),
  reportToFileAtPostEvent     = cms.untracked.string("| gzip -c > igdqm.%I.gz")
)
"""
process.source = cms.Source("PoolSource",
# replace 'myfile.root' with the source file you want to use 
fileNames = cms.untracked.vstring(
#'file:/afs/cern.ch/user/g/gpizzati/CMSSW_12_2_0_pre3/src/step2.root' 
'/store/relval/CMSSW_12_3_0_pre5/RelValTTbar_14TeV/GEN-SIM-RECO/PU_123X_mcRun4_realistic_v4_2026D88PU200-v1/10000/42b8183f-6599-4aad-b6db-f51740af9860.root'
#'/store/relval/CMSSW_12_3_0_pre5/RelValTTbar_14TeV/GEN-SIM-RECO/123X_mcRun4_realistic_v4_2026D88noPU-v1/10000/7cbdf153-c397-410e-8af5-dc3905565ced.root'
)
)
"""
process.ThroughputService = cms.Service('ThroughputService',
    eventRange = cms.untracked.uint32(10),
    eventResolution = cms.untracked.uint32(1),
    printEventSummary = cms.untracked.bool(True),
    enableDQM = cms.untracked.bool(True),
    dqmPathByProcesses = cms.untracked.bool(False),
    dqmPath = cms.untracked.string('Throughput'),
    timeRange = cms.untracked.double(1000),
    timeResolution = cms.untracked.double(1)
)

process.MessageLogger.cerr.ThroughputService = cms.untracked.PSet(
    limit = cms.untracked.int32(10000000)
)
"""
"""
process.FastTimerService.writeJSONSummary = cms.untracked.bool(True)
process.FastTimerService.jsonFileName = cms.untracked.string('resources.json')
"""
process.demo = offlinePrimaryVertices
"""
process.demo = cms.EDProducer(
    'gpuTest',
    TrackLabel = cms.InputTag("generalTracks"),
    beamSpotLabel = cms.InputTag("offlineBeamSpot"),

)
"""

process.out = cms.OutputModule("PoolOutputModule", 
    fileName= cms.untracked.string("file:/eos/user/g/gpizzati/prim-vertex/test3.root")
)

process.p = cms.Path(process.demo)
process.ep = cms.EndPath(process.out)
