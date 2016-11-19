import os
channel=["ee","emu","mumu"]

#Launch the full production
#process_bkg = ["zjjBig","zjjhpT","zjj","wwjj","zw","ttbarlep","zz","twlep","ttw","ttz","ww","signal1.5PSI2.5LQemu","signal1.6PSI2.5LQemu","signal1.7PSI2.5LQemu","signal1.8PSI2.5LQemu","signal1.9PSI2.5LQemu","signal2.0PSI2.5LQemu","signal2.1PSI2.5LQemu","signal2.2PSI2.5LQemu","signal2.3PSI2.5LQemu","signal2.4PSI2.5LQemu","signal2.0PSI1.5LQemu","signal2.5PSI2.6LQemu","signal1.5PSI2.5LQee","signal2.0PSI2.5LQee","signal1.5PSI2.5LQmumu","signal2.0PSI2.5LQmumu"]
#path_bkg = [1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1]
#nSplit = [120,30,30,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
#isSplit = [1,1,1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1]
#isMerged = [1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]

#Launch the production on a subset of samples
process_bkg = ["zjjBig"]
path_bkg = [1]
nSplit = [120]
isSplit = [1]
isMerged = [1]

#os.system('mkdir -p submitDir/')
#os.system('cd submitDir/')
for i in range(0,len(process_bkg)):
    print " This is the process number %d" % i
    process = process_bkg[i]
    print process
    execute_str = "root.exe fillHists.C\(\"backg\",\"pathFlag\",\"channel\",isMerged,RunN\)"
    doMerge = isMerged[i]
    path = 'path'+str(path_bkg[i])
    isSplited = isSplit[i]
    for j in range(0,3):
        chan = channel[j]
        filename = 'submit_'+process+'_'+chan+'.sh'
        print filename
        target = open(filename, 'w+')
        target.write("#!/bin/bash\n")
        jobname = '#PBS -N analyse_'+chan+'_'+process+'\n'
        output_jobname = '#PBS -o analyse_'+chan+'_'+process+'_out\n'
        error_jobname = '#PBS -e analyse_'+chan+'_'+process+'_err\n'
        execute_str = 'root.exe fillHists.C\(\\"'+process+'\\",\\"'+path+'\\",\\"'+chan+'\\",'+str(doMerge)+','+str(isSplited)+'\)'
        skip_all = False
        if isSplited == 1:
           #os.system('mkdir -p submitDir/'+process+'/')
           #os.system('cd submitDir/'+process+'/')           
           for k in range(0,nSplit[i]):
               filename = 'submit_'+process+'_'+chan+'_'+str(k)+'.sh'
               print filename
               target = open(filename, 'w+')
               target.write("#!/bin/bash\n")
               jobname = '#PBS -N analyse_'+chan+'_'+process+'_'+str(k)+'\n'
               output_jobname = '#PBS -o analyse_'+chan+'_'+process+'_'+str(k)+'_out\n'
               error_jobname = '#PBS -e analyse_'+chan+'_'+process+'_'+str(k)+'_err\n'
               execute_str = 'root.exe fillHists.C\(\\"'+process+'\\",\\"'+path+'\\",\\"'+chan+'\\",'+str(doMerge)+','+str(k)+'\)'
               target.write(jobname)
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
               os.system('qsub '+filename+'')
        else :
           target.write(jobname)
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
           os.system('qsub '+filename+'')   
#os.system('cd ..')
