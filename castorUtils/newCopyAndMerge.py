import os

#output
#targetDir = '/data01/sguazz'
#targetDir = '/raid/sguazz/MinBiasMC_397_altTrig_DG_10MEvt_seedSingleTrack_wEcaps'
#targetDir = '/raid/sguazz/MinBiasMC_397_altTrig_DG_10MEvt_seedSingleTrack_wEcaps_inCKFandGSF'
targetDir = '/raid/sguazz'
#targetDir = '/tmp/sguazz'
#targetDir = '/tmp'
#targetDir = '/afs/cern.ch/user/s/sguazz/scratch1'

tmpCopyDir = '/tmp/sguazz'

rootName = 'ntuple_'
rootType = 'conversion_'
#rootType = 'nuclint_'
rootSuffix = ''
#rootSuffix = 'MinBiasMC_'
#rootSuffix = 'Run2010A_'
#rootSuffix = 'HIRun2010_'
rootExt = '.root'
numberOfFilesPerChunk = 10
#numberOfFilesPerChunk = 17

#castor
#castorNumberOfFiles = 270
#castorNumberOfFiles = 92
castorNumberOfFiles = 255
#castorNumberOfFiles = 17
#castorPath = '/castor/cern.ch/user/s/sguazz/MaterialNtuples/workDir/MinBiasMC_altTrig_DG_STD/'
#castorPath = '/castor/cern.ch/user/s/sguazz/MaterialNtuples/workDir/MinBiasMC_altTrig_DG_onlyInCKFSpecialTag/'
#castorPath = '/castor/cern.ch/user/s/sguazz/MaterialNtuples/workDir/MinBiasMC2T_altTrig/'
#castorPath = '/castor/cern.ch/user/s/sguazz/MaterialNtuples/workDir/MinBiasMC_Quad'
castorPath = '/castor/cern.ch/user/s/sguazz/MaterialNtuples/workDir/MinBiasMC_i2c'
#castorPath = '/castor/cern.ch/user/s/sguazz/MaterialNtuples/workDir/MinBiasMC_i2c_sLeg'
#castorPath = '/castor/cern.ch/user/s/sguazz/MaterialNtuples/workDir/HIRun2010/'
#castorPath = '/castor/cern.ch/user/s/sguazz/MaterialNtuples/workDir/MinBiasMC_altTrig_DG_onlyInGSF/'
#castorPath = '/castor/cern.ch/user/s/sguazz/MaterialNtuples/workDir/MinBiasMC_altTrig_DG_inCKFandGSF/'
#castorPath = '/castor/cern.ch/user/s/sguazz/MaterialNtuples/workDir/MinBiasMC_altTrig_DG/'
#castorPath = '/castor/cern.ch/user/s/sguazz/MaterialNtuples/workDir/MinBiasMC_altTrig_DG_tinyTest/'
#castorPath = '/castor/cern.ch/user/s/sguazz/MaterialNtuples/workDir/MinBiasMC_altTrig_noServer/'
#castorPath = '/castor/cern.ch/user/s/sguazz/MaterialNtuples/workDir/MinBiasMC_altTrig/'
#castorPath = '/castor/cern.ch/user/s/sguazz/MaterialNtuples/workDir/Run2010A_altTrig/'
castorName = rootName
castorType = rootType
#castorSuffix = rootSuffix+'altTrig_DG_'
#castorSuffix = rootSuffix+'altTrig_'
castorSuffix = rootSuffix
castorExt = '_1.root'

#####

#addCommand = 'hadd'
addCommand = 'hadd -f'
#bkg = ' '
bkg = ' &'
doCastor = 1
doViaTMPCopy = 0
dryRun = 0

####################################################
nFile = 0
nChunk = 0
while nFile<castorNumberOfFiles :
    nChunk = nChunk+1
    chunkName = str(nChunk)
    nInChunk = 0
    mergeCommand = addCommand+' '+targetDir+'/'+rootName+rootType+rootSuffix+chunkName.zfill(2)+rootExt
    while nInChunk<numberOfFilesPerChunk :
        nFile = nFile+1
        nInChunk = nInChunk+1
        if nFile<castorNumberOfFiles+1:
            castorFileName = castorName+castorType+castorSuffix+str(nFile)+castorExt
            castorFullPath = castorPath+'/'+castorFileName
            command = 'stager_qry -M '+castorFullPath
            print command
            check = os.system(command)
            if check:
                print 'Skipping non existing file: '+castorFullPath
            else:
                if doCastor:
                    command = 'stager_get -M '+castorFullPath
                    print command
                    if not dryRun:
                        os.system(command)
                if doViaTMPCopy:
                    copyFullPath = tmpCopyDir+'/'+castorFileName
                    command = 'rfcp '+castorFullPath+' '+copyFullPath
                    print command
                    copyOk = 0
                    if not dryRun:
                        copyOk = os.system(command)
                    if copyOk:
                        print ' rfcp failed; skippin in merge: '+castorFullPath
                    else:
                        mergeCommand = mergeCommand+' '+copyFullPath
                else:
                    mergeCommand = mergeCommand+' rfio:'+castorFullPath
    mergeCommand = mergeCommand+bkg
    print '--------------------------------------------------------'
    print mergeCommand
    print '--------------------------------------------------------'
    if not dryRun:
        os.system(mergeCommand)
# Fine
print 'Done.'
