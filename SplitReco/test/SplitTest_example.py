import FWCore.ParameterSet.Config as cms

process = cms.Process("SplitReco")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.load("Configuration.StandardSequences.Services_cff")
process.load("Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff")
process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("SimTracker.TrackAssociation.TrackAssociatorByHits_cfi")

process.GlobalTag.globaltag = "START44_V9B::All"

readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring()
process.source = cms.Source ('PoolSource', fileNames=readFiles, secondaryFileNames=secFiles)

readFiles.extend( [
    'file:ciao.root'
    #'file:/raid/sguazz/2012_SummerDev/QCD_Pt-15to30_TuneZ2_7TeV_pythia6_GEN-SIM-RECODEBUG_PU_S6_START44_V5-v1.root'
    #'/store/mc/Fall11/QCD_Pt-15to30_TuneZ2_7TeV_pythia6/GEN-SIM-RECODEBUG/PU_S6_START44_V5-v1/0000/18778F04-D425-E111-8F1F-485B39800BB3.root'
    #'/store/mc/Fall11/MinBias_TuneZ2_7TeV-pythia6/GEN-SIM-RECODEBUG/NoPileUp_START44_V9B-v2/0000/0237063A-B070-E111-B606-008CFA00038C.root'
    ] );


secFiles.extend( [
    ] );

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(30) )
process.MessageLogger.cerr.FwkReport.reportEvery = 1

process.splitreco = cms.EDAnalyzer('SplitReco',
                                   	                              isMC = cms.untracked.bool( True ),
                                   #isMC = cms.untracked.bool( False ),
                                   tracks = cms.untracked.InputTag('generalTracks'),
                                   #                             tp = cms.untracked.InputTag('mergedtruth:MergedTrackTruth'),
                                   tp = cms.untracked.InputTag("mergedtruth","MergedTrackTruth"),
                                   #                             myDebug = cms.untracked.bool( False ),
                                   myDebug = cms.untracked.bool( True ),
                                   superSplit = cms.untracked.bool( False ),
                                   outputFileName = cms.untracked.string("SplitReco_test.root"),
                                   splitTrackEffHits = cms.untracked.int32(4),
                                   minSplits = cms.untracked.int32(3),
                                   kFactor = cms.untracked.double(.01),
                                   #Error treatment Barrel
                                   specialRescale = cms.untracked.bool( True ),
                                   # 0 = default: q/p, lambda, phi, dxy, dsz
                                   # 1 = q/pt, lambda, phi, dxy, dsz
                                   # 2 = q/pz, lambda, phi, dxy, dsz
                                   matrixView = cms.untracked.int32( 0 ),
                                   qOverPComp = cms.untracked.bool( False ),
                                   lambdaComp = cms.untracked.bool( True ),
                                   phiComp = cms.untracked.bool( False ),
                                   dxyComp = cms.untracked.bool( False ),
                                   dszComp = cms.untracked.bool( True ),
                                   #Error treatment Endcap
                                   EndCapSpecialRescale = cms.untracked.bool( True ),
                                   # 0 = default: q/p, lambda, phi, dxy, dsz
                                   # 1 = q/pt, lambda, phi, dxy, dsz
                                   # 2 = q/pz, lambda, phi, dxy, dsz
                                   EndCapMatrixView = cms.untracked.int32( 0 ),
                                   EndCapQOverPComp = cms.untracked.bool( False ),
                                   EndCapSLambdaComp = cms.untracked.bool( True ),
                                   EndCapPhiComp = cms.untracked.bool( True ),
                                   EndCapDxyComp = cms.untracked.bool( True ),
                                   EndCapDszComp = cms.untracked.bool( True ),
                                   #                                   
                                   etaBarrelEndcapCut =  cms.untracked.double(0.9),
                                   #
                                   # ptMinCut = cms.untracked.double(0.8),
                                   # ptMaxCut = cms.untracked.double(5.2),
                                   ptMinCut = cms.untracked.double(0.9),
                                   ptMaxCut = cms.untracked.double(1.1),
                                   nHitMinCut = cms.untracked.int32(6),
                                   etaMinCut = cms.untracked.double(0.8),
                                   etaMaxCut = cms.untracked.double(1.2),
                                   #
                                   )

process.load("Configuration.EventContent.EventContent_cff")
process.RECOEventContent.outputCommands.append('drop *')
process.RECOEventContent.outputCommands.append('keep *_*_*_SplitReco')
process.RECO = cms.OutputModule("PoolOutputModule",
                                process.RECOEventContent,
                                fileName = cms.untracked.string('reco.root')
                                )

process.p = cms.Path(process.splitreco)

