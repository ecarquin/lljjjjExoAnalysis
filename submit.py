import os
channel=["ee","emu","mumu"]
#channel=["ee"]

#Launch the full production
#process_bkg = ["zjjBig","wwjj","zw","ttbarlep","zz","twlep","ttw","ttz","ww","signal1.5PSI2.5LQemu","signal1.6PSI2.5LQemu","signal1.7PSI2.5LQemu","signal1.8PSI2.5LQemu","signal1.9PSI2.5LQemu","signal2.0PSI2.5LQemu","signal2.1PSI2.5LQemu","signal2.2PSI2.5LQemu","signal2.3PSI2.5LQemu","signal2.4PSI2.5LQemu","signal2.0PSI1.5LQemu","signal2.5PSI2.6LQemu","signal1.5PSI2.5LQee","signal2.0PSI2.5LQee","signal1.5PSI2.5LQmumu","signal2.0PSI2.5LQmumu"]
#path_bkg    = [       1,     1,   0,         0,   0,      0,    0,    0,   0,                     0,                     0,                     0,                     0,                     0,                     1,                     1,                     1,                     1,                     1,                     1,                     1,                    1,                    1,                      1,                      1]
#nSplit      = [      20,    50,  40,       100,  50,     20,    5,    5, 100,                    20,                    20,                    20,                    20,                    20,                    20,                    20,                    20,                    20,                    20,                    20,                    20,                   20,                   20,                     20,                     20]
#isSplit     = [       1,     1,   1,         1,   1,      1,    1,    1,   1,                     1,                     1,                     1,                     1,                     1,                     1,                     1,                     1,                     1,                     1,                     1,                     1,                    1,                    1,                      1,                      1]
#isMerged    = [       1,     0,   0,         0,   0,      0,    0,    0,   0,                     0,                     0,                     0,                     0,                     0,                     0,                     0,                     0,                     0,                     0,                     0,                     0,                    0,                    0,                      0,                      0]
#nFile_s     = [      20,    -1,  -1,        -1,  -1,     -1,   -1,   -1,  -1,                    -1,                    -1,                    -1,                    -1,                    -1,                    -1,                    -1,                    -1,                    -1,                    -1,                    -1,                    -1,                   -1,                   -1,                     -1,                     -1]

#Launch the full production
#process_bkg = ["wwjj","zw","ttbarlep","zz","twlep","ttw","ttz","ww","signal1.5PSI2.5LQemu","signal1.6PSI2.5LQemu","signal1.7PSI2.5LQemu","signal1.8PSI2.5LQemu","signal1.9PSI2.5LQemu","signal2.0PSI2.5LQemu","signal2.1PSI2.5LQemu","signal2.2PSI2.5LQemu","signal2.3PSI2.5LQemu","signal2.4PSI2.5LQemu","signal2.0PSI1.5LQemu","signal2.5PSI2.6LQemu","signal1.5PSI2.5LQee","signal2.0PSI2.5LQee","signal1.5PSI2.5LQmumu","signal2.0PSI2.5LQmumu"]
#path_bkg = [1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1]
#nSplit = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
#isSplit = [-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1]
#isMerged = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]

#Launch the production on a subset of samples
process_bkg = ["zjjBig"]
path_bkg    = [       1]
nSplit      = [      20]
isSplit     = [       1]
isMerged    = [       1]
nFile_s     = [      20]

#Launch the production on a subset of samples
#process_bkg = ["ttbarlep"]
#path_bkg = [0]
#nSplit   = [100]
#isSplit  = [1]
#isMerged = [0]
#isMultiF = [0]
#nFile_s  = [-1]

jname_prefix = '#PBS -N '
ename_prefix = '#PBS -e '
oname_prefix = '#PBS -o '
local_path   = '/user/e/edson/private/MG5_aMC_v2_3_3/Delphes/'

#os.system('rm -rf submitDir/')
#os.system('mkdir submitDir/')

def writeFile(filename,execute_str,jobname):
    target = open(filename, 'w+')
    target.write("#!/bin/bash\n")
    jobname_final = jname_prefix+jobname+'\n'
    output_jobname = oname_prefix+jobname+'_out\n'
    error_jobname = ename_prefix+jobname+'_err\n'
    target.write(jobname_final)
    target.write(output_jobname)
    target.write(error_jobname)
    target.write("#PBS -l mem=400mb\n")
    target.write("#PBS -l cput=1000:00:00\n")
    target.write("#PBS -l walltime=1000:00:00\n")
    target.write("#PBS -m bea\n")
    target.write("#PBS -M edson.carquin@usm.cl\n")
    target.write("cd $PBS_O_WORKDIR\n")
    target.write("\n")
    target.write("use root5\n")
    target.write("use gcc49\n")
    target.write("\n")
    target.write(execute_str)
    target.close()

for i in range(0,len(process_bkg)):
    print " This is the process number %d" % i
    process      = process_bkg[i]
    print process
    doMerge      = isMerged[i]
    #doMultiFiles = isMultiF[i]
    nFiles       = nFile_s[i]
    path         = 'path'+str(path_bkg[i])
    os.chdir(local_path)
    os.system('rm -rf submitDir/submitDir_'+process+'/')
    os.system('mkdir -p submitDir/submitDir_'+process+'/')
    os.chdir(local_path+'submitDir/submitDir_'+process+'/')
    for j in range(0,len(channel)):
        chan = channel[j]
        execute_str_init  = 'root.exe -q -b fillHists.C\(\\"'+process+'\\",\\"'+path+'\\",\\"'+chan+'\\",'+str(doMerge)+','+str(nFiles)+','
        filename_init     = 'submit_'+process+'_'+chan+'_'
        jobname_init      = 'analyse_'+chan+'_'+process+'_'
        execute_str = execute_str_init
        filename    = filename_init
        jobname     = jobname_init
        if nFiles != -1:
           for l in range(0,nFiles):
               for k in range(0,nSplit[i]):
                   os.system('mkdir submitDir_'+chan+'_'+str(l)+'_'+str(k)+'/')
                   os.system('cp '+local_path+'fillHists.C submitDir_'+chan+'_'+str(l)+'_'+str(k)+'/')
                   os.chdir('submitDir_'+chan+'_'+str(l)+'_'+str(k)+'/')
                   fname_app = str(l)+'_'+str(k)+'.sh'
                   filename += fname_app
                   print filename
                   ex_app = str(l)+','+str(k)+','+str(nSplit[i])+'\)'
                   execute_str += ex_app  
                   print execute_str
                   jname_app = str(l)+'_'+str(k)
                   jobname += jname_app
                   print jobname
                   writeFile(filename,execute_str,jobname)
                   os.system('qsub '+filename+'')
                   execute_str	= execute_str_init
                   filename    = filename_init
                   jobname     = jobname_init         
                   os.chdir(local_path+'submitDir/submitDir_'+process+'/')          
        else :           
           for k in range(0,nSplit[i]):
               os.system('mkdir submitDir_'+chan+'_'+str(k)+'/')
               os.system('cp '+local_path+'fillHists.C submitDir_'+chan+'_'+str(k)+'/')
               os.chdir('submitDir_'+chan+'_'+str(k)+'/') 
               fname_app = str(k)+'.sh'
               filename += fname_app
               print filename
               ex_app = '-1,'+str(k)+','+str(nSplit[i])+'\)'
               execute_str += ex_app
       	       print execute_str
       	       jname_app = str(k)
       	       jobname += jname_app
       	       print jobname
               writeFile(filename,execute_str,jobname)
               os.system('qsub '+filename+'')
               execute_str	= execute_str_init
               filename    = filename_init
               jobname     = jobname_init
               os.chdir(local_path+'submitDir/submitDir_'+process+'/')
