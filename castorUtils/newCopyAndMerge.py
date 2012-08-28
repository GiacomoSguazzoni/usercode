import os

#output
targetDir = '/tmp/sguazz'

tmpCopyDir = '/tmp/sguazz/copy'

rootName = 'split_'
rootType = ''
#rootType = 'conversion_'
#rootType = 'nuclint_'
rootSuffix = ''
#rootSuffix = 'MinBiasMC_'
#rootSuffix = 'Run2010A_'
#rootSuffix = 'HIRun2010_'
rootExt = '.root'
numberOfFilesPerChunk = 20
#numberOfFilesPerChunk = 17

#castor
castorNumberOfFiles = 244
castorPath = '/castor/cern.ch/user/s/sguazz/MaterialNtuples/workDir/split44'
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
#
            command = 'stager_qry -M '+castorFullPath
            print command
            fin,fout = os.popen4(command)
            result = fout.read()
            print result
            if ( 'No such file or directory' in result ) & ~( 'not on disk cache' in result ):
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
