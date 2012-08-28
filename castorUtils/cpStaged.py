import os
import subprocess

#output
targetDir = '/raid/sguazz'

tmpCopyDir = '/raid/sguazz'

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
            command = 'stager_qry -M '+castorFullPath
            print command
#            check = os.system(command)
            proc = subprocess.Popen(command,shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
            check = proc.wait()
            # Read from pipes
            check = 0;
            for line in proc.stdout:
                print("stdout: " + line.rstrip())
                if ( line.find('STAGED') ) != -1:
                    check = 1
            for line in proc.stderr:
                print("stderr: " + line.rstrip())
            if check:
                os.system("date")
                copyFullPath = tmpCopyDir+'/'+castorFileName
                command = 'rfcp '+castorFullPath+' '+copyFullPath
                print command
                copyOk = 0
                if not dryRun:
                    copyOk = os.system(command)
            else:
                print 'Skipping not staged file: '+castorFullPath
# Fine
print 'Done.'
