import FWCore.ParameterSet.Config as cms

VertexAssociatorByPositionAndTracks = cms.EDProducer('VertexAssociatorByPositionAndTracksProducer',
  absZ = cms.double(0.1),
  sigmaZ = cms.double(3),
  maxRecoZ = cms.double(1000),
  absT = cms.double(-1),
  sigmaT = cms.double(-1),
  maxRecoT = cms.double(-1),
  sharedTrackFraction = cms.double(0.6),
  trackAssociation = cms.InputTag('trackingParticleRecoTrackAsssociation'),
  mightGet = cms.optional.untracked.vstring
)

# VertexAssociatorByPositionAndTracks_tracks = cms.EDProducer('VertexAssociatorByPositionAndTracksProducer',
#   absZ = cms.double(0.9),
#   sigmaZ = cms.double(10),
#   maxRecoZ = cms.double(1000),
#   absT = cms.double(-1),
#   sigmaT = cms.double(-1),
#   maxRecoT = cms.double(-1),
#   sharedTrackFraction = cms.double(0.10),
#   trackAssociation = cms.InputTag('trackingParticleRecoTrackAsssociation'),
#   mightGet = cms.optional.untracked.vstring
# )
