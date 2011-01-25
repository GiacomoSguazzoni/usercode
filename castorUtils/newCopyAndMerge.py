import os

#output
#targetDir = '/data01/sguazz'
targetDir = '/raid/sguazz'
#targetDir = '/tmp'
#targetDir = '/afs/cern.ch/user/s/sguazz/scratch1'

rootName = 'ntuple_'
#rootType = 'conversion_'
rootType = 'nuclint_'
rootSuffix = 'MinBiasMC_'
#rootSuffix = 'Run2010A_'
rootExt = '.root'
numberOfFilesPerChunk = 26

#castor
#castorNumberOfFiles = 270
#castorNumberOfFiles = 92
castorNumberOfFiles = 102
castorPath = '/castor/cern.ch/user/s/sguazz/MaterialNtuples/workDir/MinBiasMC_altTrig_DG/'
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
do = 1

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
            castorFullPath = castorPath+'/'+castorName+castorType+castorSuffix+str(nFile)+castorExt
            command = 'stager_get -M '+castorFullPath
            if do:
                os.system(command)
            print command
            mergeCommand = mergeCommand+' rfio:'+castorFullPath
    mergeCommand = mergeCommand+bkg
    if do:
        os.system(mergeCommand)
    print mergeCommand
# Fine
print 'Done.'
