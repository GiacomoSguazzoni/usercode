import subprocess

#params
defaultEOSRootPath = '/store/user/sguazz/qcd1530_nH4_specRescalePzView'
defaultTmpPath = '/tmp/sguazz'
defaultTargetPath = '/afs/cern.ch/user/s/sguazz/myWorkarea/merge/'
filesPerChunk = 10
doCopy = 1

#base commands
defaultEOSlistCommand = 'cmsLs '
defaultEOSCopyCommand = 'cmsStage '
defaultHaddCommand = 'hadd '

#pro = subprocess.Popen("/afs/cern.ch/cms/caf/scripts/cmsLs /store/user/sguazz/store/user/sguazz", stdout=PIPE, stderr=PIPE)

def getFileList():
    theCommand = defaultEOSlistCommand+' '+defaultEOSRootPath
    dirList = subprocess.Popen(["/bin/sh","-c",theCommand], stdout=subprocess.PIPE)        
    list = []
    for line in dirList.stdout.readlines():
        strippedLine = line.rstrip('\n')
        if strippedLine != '':
            print strippedLine
            filename = strippedLine.split('/')[-1]
            print filename
            list.append(filename)
    return list

def parseFileList(fileList):
    haddList = []
    icount = 0
    ichunk = 0
    haddChunkCommand = ''
    for file in fileList:
        destinationPath = defaultTmpPath+'/'+file
        theCommand = defaultEOSCopyCommand+' '+defaultEOSRootPath+'/'+file+' '+destinationPath
        if ( doCopy ):
            subprocess.call(["/bin/sh","-c",theCommand], stdout=subprocess.PIPE)        
        print theCommand
        if ( icount%filesPerChunk ):
            haddChunkCommand += ' '+destinationPath
        else:
            if ( ichunk ):
#                print '>>>>for chunk '+str(ichunk)+': '+haddChunkCommand
                haddList.append(haddChunkCommand)                
            ichunk+=1
            haddChunkCommand = defaultHaddCommand+' '+defaultTargetPath+'/chunk_'+str(ichunk)+'.root '+destinationPath
        icount+=1
    if ( ichunk ):
#        print '>>>>for chunk '+str(ichunk)+': '+haddChunkCommand
        haddList.append(haddChunkCommand)                
    return haddList

def doMerge(haddList):
    for haddChunkCommand in haddList:
        print haddChunkCommand
        subprocess.call(["/bin/sh","-c",haddChunkCommand], stdout=subprocess.PIPE)        
        
fileList = getFileList()
haddList = parseFileList(fileList)
doMerge(haddList)


