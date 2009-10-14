import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
#process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(10000)
#process.MessageLogger.categories     = ["TrackAssociator","TrackValidator", "TestHits","TrackingRegressionTest",
#       "TrackFitters","TrackProducer","FwkSummary","RFIOFile","FwkReport","DDLParser"]
#process.MessageLogger.debugModules   = [ "KFTrajectoryFitterESProducer" ]
#process.MessageLogger.categories     = ["TrackFitters"]
#process.MessageLogger.destinations = ['log']

process.load("Configuration.StandardSequences.Services_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = "MC_31X_V3::All"
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("SimTracker.TrackAssociation.TrackAssociatorByHits_cfi")


process.source = cms.Source ("PoolSource",
                             fileNames = cms.untracked.vstring(
# 01
'/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0019/FE72E489-9190-DE11-923F-0018F3D096DC.root',
'/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0019/F81E5C6E-6F8D-DE11-B147-001731AF692D.root',
'/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0019/F0816689-9190-DE11-A00F-0018F34D0D62.root',
'/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0019/E630BAA4-9190-DE11-AADC-001BFCDBD160.root',
'/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0019/E4F97400-A38D-DE11-9E0A-001BFCDBD160.root',
'/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0019/E4E58ABC-9190-DE11-AEC3-001A92971B64.root',
'/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0019/E26CABAD-9190-DE11-A8B5-001731AF6BC1.root',
'/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0019/E08DFF12-A38D-DE11-B346-0018F3D0960E.root',
'/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0019/DA92B990-9190-DE11-AFC0-001A92810AE0.root',
'/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0019/D6848CB6-9190-DE11-AABC-0018F3D09688.root',
'/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0019/D4A5A100-A38D-DE11-AFA0-001A92810AA4.root',
'/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0019/D24300BF-9190-DE11-9AB2-001731AF6785.root',
'/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0019/D0A748BF-9190-DE11-B9D1-003048678B72.root',
'/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0019/D08BC448-6F8D-DE11-9949-00304867915A.root',
'/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0019/CEE3178B-6F8D-DE11-8757-001731AF698D.root',
'/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0019/CC9A9515-A38D-DE11-BC5B-0018F3D0970C.root',
'/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0019/CA5B22FC-A28D-DE11-A160-001BFCDBD11E.root',
'/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0019/C80E9113-A38D-DE11-B0E2-0018F3D096BC.root',
'/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0019/C4AA8DFC-A28D-DE11-BE9D-0018F3D09708.root',
'/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0019/C06E2E05-A38D-DE11-9F6A-0018F3D09608.root',
'/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0019/BEFCC14B-6F8D-DE11-B05F-001731AF66BF.root',
'/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0019/BEA587AD-9190-DE11-801F-001731AF6BCD.root',
'/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0019/B82E5B8B-6F8D-DE11-A093-001731A28A31.root',
'/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0019/B6D91F5A-6F8D-DE11-840B-0030486792BA.root',
'/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0019/B64D5810-A38D-DE11-89F5-0018F3D096F8.root',
'/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0019/B23594B8-9190-DE11-BB18-001A92971AA8.root',
'/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0019/9C85A1A3-9190-DE11-9C0D-001A92971B5E.root',
'/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0019/9C824BC5-9190-DE11-ADAB-0018F3D096DC.root',
'/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0019/9687ACF0-0B90-DE11-928A-00304867920C.root',
'/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0019/961DABFE-A28D-DE11-A43E-001A92810AA4.root',
'/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0019/945B6085-9190-DE11-A6B4-0018F3D09704.root',
'/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0019/92842085-9190-DE11-AEC8-001A92810AE4.root',
'/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0019/90E4A303-A38D-DE11-9DCD-001A92810AA8.root',
'/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0019/8231CF98-9190-DE11-8F11-0018F3D09688.root',
'/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0019/7489BAFA-A28D-DE11-83BE-0018F3D096EE.root',
'/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0019/70B0E8A3-9190-DE11-8A86-001A92971B7E.root',
'/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0019/707AD44F-6F8D-DE11-BA1F-001731AF699D.root',
'/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0019/6C888704-A38D-DE11-83E0-001A92971B0E.root',
'/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0019/6C559C9F-9190-DE11-9EC1-003048678FD6.root',
'/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0019/68528E1B-8A90-DE11-BB05-0018F3D09684.root',
'/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0019/667DEAF1-0B90-DE11-B65E-0030486792BA.root',
'/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0019/626C598C-6F8D-DE11-BDD6-001731A28998.root',
'/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0019/6209CA83-6F8D-DE11-815B-001731AF66F7.root',
'/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0019/5E28FD01-A38D-DE11-B123-001A92971B68.root',
'/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0019/5E003565-9190-DE11-B5BE-0018F3D09620.root',
'/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0019/5028DDC5-9190-DE11-B3F3-001A92971B16.root',
'/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0019/4CB531AC-9190-DE11-AEB3-001731AF65EB.root',
'/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0019/4AE86493-9190-DE11-9F70-0018F3D096DA.root',
'/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0019/4A8D0913-A38D-DE11-BCEA-001A928116B4.root',
'/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0019/48A1AC01-728E-DE11-8152-001731EF61B4.root',
'/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0019/38FE5B02-A38D-DE11-AB04-001BFCDBD19E.root',
'/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0019/34E20389-9190-DE11-B470-001BFCDBD100.root',
'/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0019/3294C94C-6F8D-DE11-ACC7-001731A28F19.root',
'/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0019/2EE615FA-0B90-DE11-B6D1-001731AF66F5.root',
'/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0019/2E805D26-6F8D-DE11-B261-0018F3D09630.root',
'/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0019/28F2416A-9190-DE11-BB0D-001731AF6A4B.root',
'/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0019/28A8CC5D-6F8D-DE11-8E98-001A92810ADE.root',
'/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0019/26DC9F15-A38D-DE11-9835-0018F3D095FA.root',
'/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0019/1CEFAD03-A38D-DE11-B5A8-0018F3C3E3A6.root',
'/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0019/169E438D-9190-DE11-A1CA-0018F3D09700.root'
,
# 02
'/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0019/0682C048-6F8D-DE11-A83E-00304867915A.root',
'/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0019/00094F00-A38D-DE11-8562-001BFCDBD19E.root',
'/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0018/FC79166C-538D-DE11-9B03-001A92971B3C.root',
'/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0018/F81FC282-538D-DE11-B178-0018F3D09616.root',
'/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0018/F2625E9A-538D-DE11-9E70-001A92971BB8.root',
'/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0018/EC7A555A-538D-DE11-9AF6-0018F3D09660.root',
'/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0018/D4BCF347-528D-DE11-8B50-003048679296.root',
'/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0018/C0C9B13B-528D-DE11-8423-001A92971B94.root',
'/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0018/BE22CF52-538D-DE11-AA9D-0018F3D095EA.root',
'/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0018/BC426ED3-518D-DE11-BAC0-0018F3D09654.root',
'/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0018/A233EA88-538D-DE11-9DAB-0018F3D09628.root',
'/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0018/A0F65C96-538D-DE11-A1BD-001731AF687F.root',
'/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0018/A030AE6A-538D-DE11-8297-0018F3D09616.root',
'/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0018/9A3EA454-0C8D-DE11-BCE8-0018F3D096BA.root',
'/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0018/94056E80-648D-DE11-A808-0018F3D095EC.root',
'/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0018/88FEF25C-0C8D-DE11-9E2D-0018F3D096CA.root',
'/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0018/86179BA9-5F8D-DE11-8851-0018F3D09688.root',
'/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0018/84352522-508D-DE11-B032-001731AF68AB.root',
'/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0018/7C758167-538D-DE11-9171-0018F3D09704.root',
'/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0018/781FC744-538D-DE11-9247-001731EF61B4.root',
'/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0018/766F37F1-518D-DE11-9F44-0018F3D096A6.root',
'/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0018/64265E8D-538D-DE11-BD01-0018F3D0965E.root',
'/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0018/600B034F-528D-DE11-B7E8-0018F3D096E0.root',
'/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0018/56FF22DC-518D-DE11-8448-0018F3D0966C.root',
'/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0018/48B4AE7A-538D-DE11-B8B2-001731AF66F7.root',
'/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0018/3A44AB8E-0C8D-DE11-BA6E-0018F3D095EA.root',
'/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0018/38753A57-538D-DE11-B8C9-0018F3D0965E.root',
'/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0018/2CA0EF51-0C8D-DE11-9FF4-0018F3D09704.root',
'/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0018/2425BB77-0C8D-DE11-974E-0018F3D096A2.root',
'/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0018/0E1D3877-538D-DE11-A539-001A92971B9A.root',
'/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0018/027306CF-0B8D-DE11-80B6-001731AF684B.root',
'/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0017/E086A9E4-FC8C-DE11-8843-00304867906C.root',
'/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0017/CE7978C3-F58C-DE11-B3F4-001731AF66A1.root',
'/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0017/C0561911-068D-DE11-9C61-001A92971BB8.root',
'/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0017/AC23AEFD-F88C-DE11-99E5-0018F3D09636.root',
'/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0017/9088DE7C-058D-DE11-908E-001731AF65E9.root',
'/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0017/903DEA1A-068D-DE11-8B90-0018F3D09704.root',
'/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0017/7CAB3520-FC8C-DE11-90F5-001BFCDBD154.root',
'/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0017/7C01B2C4-FD8C-DE11-8DFE-001A92811716.root',
'/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0017/6CCDE97B-028D-DE11-8E62-001731AF66A1.root',
'/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0017/66457ACF-F78C-DE11-B84C-0018F3D095FC.root',
'/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0017/561B7BE9-028D-DE11-BB49-003048678FC4.root',
'/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0017/52C61B56-FE8C-DE11-95F4-0018F3D0969A.root',
'/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0017/52016110-FB8C-DE11-95AD-0018F3D0968C.root',
'/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0017/4E45B5F8-0A8D-DE11-93A8-001731A281B1.root',
'/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0017/46C0D095-FD8C-DE11-BDF5-001731AF684B.root',
'/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0017/429A2B65-028D-DE11-966A-0017312B577F.root',
'/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0017/261C4D07-FC8C-DE11-92E5-001A92810AE4.root',
'/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0017/1E891BF7-0A8D-DE11-9D96-003048679228.root',
'/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0017/1CEA3CE3-FD8C-DE11-8577-0018F3D09696.root',
'/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0017/186E54B6-FC8C-DE11-B808-0030486790A6.root',
'/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0017/1660BAC9-F78C-DE11-802F-001A9281170E.root',
'/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0017/0E3AA430-FD8C-DE11-9C42-003048678ADA.root',
'/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0017/0A3BF776-0C8D-DE11-A484-0018F3D095FE.root',
'/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0017/0849F9CE-0A8D-DE11-87C9-00304867915A.root',
'/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0015/F4E8FE60-928C-DE11-A704-0018F3D09630.root',
'/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0015/EAD55E64-928C-DE11-A38F-001731AF68CF.root',
'/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0015/E2AE1465-928C-DE11-8283-001731A28BE1.root',
'/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0015/DCAD4E73-928C-DE11-9520-001731A28A31.root',
'/store/data/CRAFT09/Cosmics/RAW-RECO/GR09_31X_V5P_SuperPointing-332_v4/0015/D8A26728-928C-DE11-85E7-00304867920C.root'
    ),
                             )


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.demo = cms.EDAnalyzer('SplitTest',
                              isMC = cms.untracked.bool( False ),
                              tracks = cms.untracked.InputTag('generalTracks'),
                              cosmics = cms.untracked.InputTag('ctfWithMaterialTracksP5'),
                              simtracks = cms.untracked.InputTag('g4SimHits'),
                              tp = cms.untracked.InputTag('mergedtruth:MergedTrackTruth'),
                              myDebug = cms.untracked.bool( False ),
#                              myDebug = cms.untracked.bool( True ),
                              outputFileName = cms.untracked.string("/raid/sguazz/superPoint_part01.root"),
                              splitTrackEffHits = cms.untracked.int32(4),
                              minSplits = cms.untracked.int32(3),
                              kFactor = cms.untracked.double(.01),
                              poolCut = cms.untracked.double(7.),
                              goodruns = cms.untracked.vuint32(109035,109468,109470,109472,109474,109490,
                                                               109504,109508,109519,109524,109562,109573,
                                                               109575,109578,109586,109606,110388,110395,
                                                               110397,110409,110419,110428,110431,110437,
                                                               110440,110452,110476,110485,110496,110508,
                                                               110520,110535,110546,110683,110686,110689,
                                                               110708,110784,110832,110835,110842,110846,
                                                               110878,110894,110900,110916,110921,110924,
                                                               110958,110987,110998,111017,111023,111045,
                                                               111047,111123,111125)
                              )

process.load("Configuration.EventContent.EventContent_cff")
process.RECOEventContent.outputCommands.append('drop *')
process.RECOEventContent.outputCommands.append('keep *_*_*_Demo')
process.RECO = cms.OutputModule("PoolOutputModule",
    process.RECOEventContent,
    fileName = cms.untracked.string('reco.root')
)

process.p = cms.Path(process.demo)

#process.outpath = cms.EndPath(process.RECO)
