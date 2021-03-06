import os
channel=["ee","emu","mumu"]

#Launch the full production
process_bkg = ["zjjBig","wwjj","zw","ttbarlep","zz","twlep","ttw","ttz","ww","signal1.5PSI2.5LQemu","signal1.6PSI2.5LQemu","signal1.7PSI2.5LQemu","signal1.8PSI2.5LQemu","signal1.9PSI2.5LQemu","signal2.0PSI2.5LQemu","signal2.1PSI2.5LQemu","signal2.2PSI2.5LQemu","signal2.3PSI2.5LQemu","signal2.4PSI2.5LQemu","signal2.0PSI1.5LQemu","signal2.5PSI2.6LQemu","signal1.5PSI2.5LQee","signal2.0PSI2.5LQee","signal1.5PSI2.5LQmumu","signal2.0PSI2.5LQmumu"]
nSplit      = [      20,    50,  40,	   100,  50,     20,    5,    5, 100,                    20,                    20,                    20,                    20,                    20,                    20,                    20,                    20,                    20,                    20,                    20,                    20,                   20,                   20,                     20,                     20]
nFile_s     = [      20,    -1,  -1,        -1,  -1,     -1,   -1,   -1,  -1,                    -1,                    -1,                    -1,                    -1,                    -1,                    -1,                    -1,                    -1,                    -1,                    -1,                    -1,                    -1,                   -1,                   -1,                     -1,                     -1]

#Launch the production on a subset of samples
#process_bkg = ["ttbarlep"]
#nSplit = [100]
#nFile_s = [-1]

#Launch the production on a subset of samples
#process_bkg = ["zjjBig"]
#nSplit = [20]
#nFile_s = [20]

tag = '30_03_17/'
#os.system('mkdir '+tag+' ')
root_path = '/user/e/edson/private/MG5_aMC_v2_3_3/Delphes/'
#os.chdir('09_03_17')
#os.system('mkdir -p submitDir/')
#os.system('cd submitDir/')
for i in range(0,len(process_bkg)):
    print " This is the process number %d" % i
    process = process_bkg[i]
    print process
    proc_dir = 'submitDir/submitDir_'+process+'/'           
    os.chdir(proc_dir)
    for j in range(0,len(channel)):
        chan = channel[j]
        chan_dir = 'submitDir_'+chan+'_'
        #hadd_str_base = 'hadd -f '
        #hadd_str_out  = 'hists_'+chan+'_'+process+'.root ' 
        #hadd_str_k_full = ''
        if nFile_s[i] != -1:
            for k in range(0,nFile_s[i]):
                #hadd_str_k = 'hists_'+chan+'_'+process+'_'+str(k)+'.root '
                #path_str_l = ''
                for l in range(0,nSplit[i]):
                   kl_dir = str(k)+'_'+str(l)+'/'
                   file_name = 'hists_'+chan+'_'+process+'_'+str(k)+'_'+str(l)+'.root'
                   path_str_l = chan_dir+kl_dir+tag+file_name 
                   #print 'This is the file: ', path_str_l
                   if os.path.exists(path_str_l) == False: 
                      print 'file: ', path_str_l, ' doesnt exist, going to retry !'
                      #print chan_dir+kl_dir
                      os.chdir(chan_dir+kl_dir)
                      os.system('pwd')
                      os.system('qsub submit_'+process+'_'+chan+'_'+str(k)+'_'+str(l)+'.sh')
                   os.chdir(root_path+proc_dir)
        else:
            for k in range(0,nSplit[i]):
                k_dir = str(k)+'/'
                file_name = 'hists_'+chan+'_'+process+'_'+str(k)+'.root'
                path_str_k = chan_dir+k_dir+tag+file_name
                #print 'This is the file: ', path_str_l
                if os.path.exists(path_str_k) == False:
                   print 'file: ', path_str_k, ' doesnt exist, going to retry !'
                   #print chan_dir+kl_dir
                   os.chdir(chan_dir+k_dir)
                   os.system('pwd')
                   os.system('qsub submit_'+process+'_'+chan+'_'+str(k)+'.sh')
                os.chdir(root_path+proc_dir)
        #os.system(hadd_string)
        #os.system('rm '+'hists_'+chan+'_'+process+'_*.root')
    os.chdir(root_path)
os.system('pwd')
