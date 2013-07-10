import FWCore.ParameterSet.Config as cms

process = cms.Process("SplitReco")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.load("Configuration.StandardSequences.Services_cff")
process.load("Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff")
process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("SimTracker.TrackAssociation.TrackAssociatorByHits_cfi")

process.GlobalTag.globaltag = "START53_V15::All"

readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring()
process.source = cms.Source ('PoolSource', fileNames=readFiles, secondaryFileNames=secFiles)

readFiles.extend( [
    'file:/raid/sguazz/test/mc/TwentyMuons_Pt0d9to1d1_Summer12_DR53X-NominalGeo_NoPileUp_START53_V15-v1_GEN-SIM-RECODEBUG.root'
    #'file:/raid/sguazz/2012_SummerDev/QCD_Pt-15to30_TuneZ2_7TeV_pythia6_GEN-SIM-RECODEBUG_PU_S6_START44_V5-v1.root'
    #'/store/mc/Fall11/QCD_Pt-15to30_TuneZ2_7TeV_pythia6/GEN-SIM-RECODEBUG/PU_S6_START44_V5-v1/0000/18778F04-D425-E111-8F1F-485B39800BB3.root'
    #'/store/mc/Fall11/MinBias_TuneZ2_7TeV-pythia6/GEN-SIM-RECODEBUG/NoPileUp_START44_V9B-v2/0000/0237063A-B070-E111-B606-008CFA00038C.root'
    ] );


secFiles.extend( [
    ] );

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.ectest = cms.EDAnalyzer('ecReco',
                                isMC = cms.untracked.bool( True ),
                                #isMC = cms.untracked.bool( False ),
                                tracks = cms.untracked.InputTag('generalTracks'),
                                #                             tp = cms.untracked.InputTag('mergedtruth:MergedTrackTruth'),
                                tp = cms.untracked.InputTag("mergedtruth","MergedTrackTruth"),
                                myDebug = cms.untracked.bool( False ),
                                #myDebug = cms.untracked.bool( True ),
                                outputFileName = cms.untracked.string("ecTest.root"),
                                # ptMinCut = cms.untracked.double(0.8),
                                # ptMaxCut = cms.untracked.double(5.2),
                                ptMinCut = cms.untracked.double(0.2),
                                ptMaxCut = cms.untracked.double(10.0),
                                etaMinCut = cms.untracked.double(0.0),
                                etaMaxCut = cms.untracked.double(2.5),
                                #
                                nHitMinCut = cms.untracked.int32(8),
                                nHitTlMinCut = cms.untracked.int32(4),
                                #
                                nHitTibMinCut = cms.untracked.int32(0),
                                nHitSteTibMinCut = cms.untracked.int32(0),
                                useTIB = cms.untracked.bool( False ),
                                #
                                nHitTidMinCut = cms.untracked.int32(0),
                                nHitSteTidMinCut = cms.untracked.int32(0),
                                useTID = cms.untracked.bool( True ),
                                #
                                nHitTobMinCut = cms.untracked.int32(0),
                                nHitSteTobMinCut = cms.untracked.int32(0),
                                useTOB = cms.untracked.bool( False ),
                                #
                                nHitTecMinCut = cms.untracked.int32(0),
                                nHitSteTecMinCut = cms.untracked.int32(0),
                                useTEC = cms.untracked.bool( True ),
                                #
                                useMono = cms.untracked.bool( False ),
                                useStereo = cms.untracked.bool( False ),
                                useStereoIfGlued = cms.untracked.bool( True ),
                                #
                                # Activate to have a tracklet == full track; only for debugging
                                useAll = cms.untracked.bool( False )
#                                Builder = cms.untracked.string('WithTrackAngle')
)

process.load("Configuration.EventContent.EventContent_cff")
process.RECOEventContent.outputCommands.append('drop *')
process.RECOEventContent.outputCommands.append('keep *_*_*_ecTest')
process.RECO = cms.OutputModule("PoolOutputModule",
                                process.RECOEventContent,
                                fileName = cms.untracked.string('reco.root')
                                )

process.p = cms.Path(process.ectest)

