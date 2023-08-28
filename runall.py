import subprocess
import sys
import os
# opening the file in read mode
my_file = open(sys.argv[1], "r")
  
# reading the file
data = my_file.read()
  
# replacing end splitting the text 
# when newline ('\n') is seen.
filelist = data.split("\n")
#itdtseed = [["randitdt", "randitdt", "randseed"], ["0.5", "0.3", "maxdegree"], ["0.1", "0.3", "maxdegree"]]
itdtseed = [["0.5", "0.3", "maxdegree"], ["0.1", "0.3", "maxdegree"]]
i = 0
while(i < len(filelist)-1):   
    for idtype in itdtseed:        
        SRsize = ["3", "5", "10"]
        for rsc in SRsize:
            TSsize = ["3", "5", "10"]
            for tsc in TSsize:                                                                
                noiter = 1
                #if(idtype[0] == "randitdt"):
                #    noiter = 1                    
                runcmd = "./a.out "+ filelist[i]+" "+rsc+" "+tsc+" "+idtype[0]+" "+idtype[1]+" "+idtype[2]+" runall "
                for f in range(2):
                    runcmd =  runcmd + filelist[f+1]+" "
                runcmd = runcmd + "10 50 100 2000"
                #for k in range(noiter):                        
                #print(runcmd)
                os.system(runcmd)                                                          
    i = i + 3     
