import os

#output
#targetDir = '/data01/sguazz'
targetDir = '/raid/sguazz'
#targetDir = '/tmp/sguazz'
#targetDir = '/tmp'
#targetDir = '/afs/cern.ch/user/s/sguazz/scratch1'

tmpCopyDir = '/tmp/sguazz'

rootName = 'ntuple_'
rootType = 'conversion_'
#rootType = 'nuclint_'
rootSuffix = 'MinBiasMC_'
#rootSuffix = 'Run2010A_'
rootExt = '.root'
numberOfFilesPerChunk = 26

#castor
#castorNumberOfFiles = 270
#castorNumberOfFiles = 92
castorNumberOfFiles = 102
castorPath = '/castor/cern.ch/user/s/sguazz/MaterialNtuples/workDir/MinBiasMC_altTrig_DG_STD/'
#castorPath = '/castor/cern.ch/user/s/sguazz/MaterialNtuples/workDir/MinBiasMC_altTrig_DG_onlyInGSF2/'
#castorPath = '/castor/cern.ch/user/s/sguazz/MaterialNtuples/workDir/MinBiasMC_altTrig_DG_inCKFandGSF/'
###castorPath = '/castor/cern.ch/user/s/sguazz/MaterialNtuples/workDir/MinBiasMC_altTrig_DG/'
#castorPath = '/castor/cern.ch/user/s/sguazz/MaterialNtuples/workDir/MinBiasMC_altTrig_DG_tinyTest/'
#castorPath = '/castor/cern.ch/user/s/sguazz/MaterialNtuples/workDir/MinBiasMC_altTrig_noServer/'
#castorPath = '/castor/cern.ch/user/s/sguazz/MaterialNtuples/workDir/MinBiasMC_altTrig/'
#castorPath = '/castor/cern.ch/user/s/sguazz/MaterialNtuples/workDir/Run2010A_altTrig/'
castorName = rootName
castorType = rootType
castorSuffix = rootSuffix+'altTrig_DG_'
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
