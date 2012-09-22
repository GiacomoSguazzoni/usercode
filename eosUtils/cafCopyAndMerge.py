import os

#output
targetDir = '/raid/sguazz'
tmpCopyDir = '/tmp/sguazz'
tmpExt = '.root'

chunkName = 'ntuple_nuclint_'
chunkExt = '.root'
numberOfFilesPerChunk = 10

#eos
eosNumberOfFiles = 86
eosPath = '/store/caf/user/sguazz/Run2011A_2T_162713_MB1_rereco/'


#####

addCommand = 'hadd -f'
bkg = ' &'
dryRun = 0


# Empty tmp dir (to be on the safe side)
command = 'rm -Rf '+tmpCopyDir+'/*'
print command
if not dryRun:
    os.system(command)

####################################################
nFile = 0
nChunk = 0
# while nFile<eosNumberOfFiles :

#
# Get list of files in dir
#
command = "cmsLs "+eosPath+" | grep root | awk '{ print $5 }'"
#fileList = os.system(command)

nInChunk = 0
mergeCommand = addCommand+' '+targetDir+'/'+chunkName+str(nChunk).zfill(2)+chunkExt

for myFile in os.popen(command).readlines() :
    eosFullPath=myFile.rstrip()
##  while fileList :
    nInChunk = nInChunk+1
    nFile = nFile+1

    copyFullPath = tmpCopyDir+'/sguazzFile_'+str(nFile)+tmpExt
    command = 'cmsStage '+eosFullPath+' '+copyFullPath
    print command
    if not dryRun:
        os.system(command)

    mergeCommand = mergeCommand + ' ' + copyFullPath

    if nInChunk==numberOfFilesPerChunk :
        mergeCommand = mergeCommand+bkg
        print mergeCommand
        if not dryRun:
            os.system(mergeCommand)

        nInChunk = 0 
        nChunk = nChunk+1
        mergeCommand = addCommand+' '+targetDir+'/'+chunkName+str(nChunk).zfill(2)+chunkExt

# Flush last chunk
mergeCommand = mergeCommand+bkg
print mergeCommand
if not dryRun:
    os.system(mergeCommand)


# Fine
print 'Done.'
