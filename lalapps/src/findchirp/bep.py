#!/usr/bin/python2.2 
"""
python script to create a condor_script and xml prototype  
arguments might be changed directly in this file.
"""

__author__ = 'Thomas Cokelaer <Thomas.Cokelaer@astro.cf.ac.uk>'
__date__ = '$Date$'
__version__ = '$Revision$'


import sys
import os
import optparse
import locale
import math
from optparse import OptionParser
import time  
import ConfigParser

#find the path for lalapps_bankefficiency
user = os.getlogin()
uname = os.uname()
host = uname[1]
#if host.find('explorer')>=0:



def set_predefined_search_parameter(BE):
	if BE['search']=='BNS':
	    BE['bank-mass-range'] = '1 3'
	    BE['signal-mass-range'] = '1 3'
	    BE['sampling'] = 4096
	    BE['bank-ffinal']= BE['sampling']/2 - 1
    
	elif BE['search']=='BBH':
	    BE['bank-mass-range'] = '3 30'
	    BE['signal-mass-range'] = '3 30'
	    BE['sampling'] = 4096
	    BE['bank-ffinal'] = BE['sampling']/2 - 1
	
	elif BE['search']=='BHNS':
	    BE['bank-mass-range'] = '1 63'
	    BE['signal-mass-range'] = '1 60'
	    BE['sampling'] = 4096
	    BE['bank-ffinal'] = BE['sampling']/2 - 1
	
	elif BE['search']=='S5':
	    BE['bank-mass-range'] = '1 30'
	    BE['signal-mass-range'] = '1 30'
	    BE['sampling'] = 4096
	    BE['bank-ffinal'] = BE['sampling']/2 - 1
	
	elif BE['search']=='BCV':
	    BE['bank-psi0-range'] = '10000 550000'
	    BE['bank-psi3-range'] = '-5000 -10'
	    BE['bank-number-fcut'] = '3'
	    BE['signal-mass-range'] = '3 80'
	    BE['bank-mass-range'] = '3 80'
	    BE['signal-max-total-mass'] =  80
	    BE['bank-alpha'] = 0.01
	    BE['sampling'] = 4096
	    BE['bank-ffinal'] = BE['sampling']/2 - 1
	    BE['template'] = 'BCV'
	    BE['bank-inside-polygon'] = 1
	   
	   
	elif BE['search']=='PBH':
	    BE['bank-mass-range'] = '.3 1'
	    BE['signal-mass-range'] = '.3 1'
	    BE['sampling'] = 4096
	    BE['bank-ffinal'] = BE['sampling']/2 - 1
	else:
	    BE['bank-mass-range'] = '1 3'
	    BE['signal-mass-range'] = '1 3'
	    BE['sampling'] = 4096
	    BE['bank-ffinal'] = BE['sampling']/2 - 1
	return BE

def create_condor_file(configcp):
  """
  """
  fp = open('bep.sub','w');
  fp.write('Executable   = ' +configcp.get("main", "executable")+"\n")
  fp.write('Universe     = vanilla\n')

  if host.find('coma')>=0:        
    fp.write('Environment  = LD_LIBRARY_PATH=/usr/lscsoft/non-lsc/lib\n')
     
    #fp.write('Requirements = Memory >=128 && OpSys == "LINUX" && FileSystemDomain == "explorer" && UidDomain == "explorer"\n')
    #fp.write('+MaxHours =40\n\n')

  arguments = ""
  for i in configcp.items("general"):
    arguments = arguments + ' --'+i[0] +' ' +i[1]
  for i in configcp.items("bank"):
    arguments = arguments + ' --'+i[0] +' ' +i[1]
  for i in configcp.items("signal"):
    arguments = arguments + ' --'+i[0] +' ' +i[1]
    
  n = float(configcp.get("simulation", "ntrial"))
  N = float(configcp.get("simulation", "njobs"))
  print n/N
  
  
  arguments = arguments + ' --ntrial '+str( int(n/N))



  fp.write('Arguments = ' + arguments + ' --seed $(macroseed)')
  fp.write('\n\n')
  fp.write('priority = 10\n')
  
  
  if host.find('ldas-grid')>=0:        
    index = 1
    this_path = '/usr1/'+user+'/'
    this_file = 'tmp'+str(index)
    while os.path.isfile(this_path+this_file)==True:
      index=index+1
      this_file = 'tmp'+str(index)
			
    msg = 'log = '+this_path+this_file+'\n'	
  else:
    msg = 'log = ./log/tmp\n'
    
  fp.write(msg)
  msg = 'output = ./out_$(macroseed) \n'
  fp.write(msg)
  msg = 'error = ./log/err_$(macroseed) \n'
  fp.write(msg)
  fp.write('notification = never\n')

  fp.write('Queue 1')
  fp.close()
  
  return arguments

def create_bank(configcp, arguments):
  """
  
  """
  arguments = arguments + ' --n 1 --check --print-bank --print-xml'
  os.system('rm -f BE_Bank.dat BE_Bank.xml')
  print '###'
  print ' We are creating the template bank for sanity check. Please wait'
  fp =open('BankEfficiency_createbank','w');
  fp.write( configcp.get("main", "executable") + arguments +' 1> bep_bank.out 2>bep_bank.err'+'\n')
  fp.close()
  os.system('chmod 755 BankEfficiency_createbank')
  a=os.system('./BankEfficiency_createbank')
  
  if a==0:
    print '... done (your parameters seems correct). See BE_Bank.xml file.'
  else:
    print '... failed (your parameters seems correct)'
    sys.exit()

def create_dag_file(configcp):
  """ 
  """
  njobs = int(configcp.get("simulation", "njobs"))
  print '--- Generating the dag file'
  fp=open('bep.dag', 'w')
  for id in range(1,njobs+1,1):
    fp.write('JOB '+str(id)+' bep.sub\n')
    fp.write('VARS '+str(id)+' macroseed="'+str(id)+'"\n')

  fp.write('JOB '+str(njobs+1)+ ' finalise.sub\n' )
  
  for id in range(1,njobs+1,1):
    fp.write('PARENT ' + str(id)+' CHILD '+str(njobs+1)+'\n')
    

                
  fp.close()
  print '... done'

def create_finalise_script(configcp):
  """
  """
  fl = configcp.get("general", "fl")
  noise_model = configcp.get("general", "noise-model")
  grid = configcp.get("bank", "bank-grid-spacing")
  mm = configcp.get("bank", "mm")
  template = configcp.get("bank", "template")
  signal = configcp.get("signal", "signal")
  
  fp = open('finalise.sh', 'w')
  fp.write('#!/bin/sh\n')
  fp.write('cp TMPLTBANK.xml BE_Bank.xml\n')
  fp.write('rm -f Trigger.dat ; find . -name "out_*" | awk \'{print "cat  " $1 ">> Trigger.dat"}\' > script.sh; chmod 755 script.sh ; ./script.sh; \n')
  fp.write(configcp.get("main", "executable") +' --ascii2xml \n')
  fp.write('mv Trigger.xml Trigger_' + noise_model +'_'+fl+'_'+grid+'_'+template+'_'+signal+'_'+mm+'.xml')
  fp.close()
  os.system('chmod 755 finalise.sh')
        

def check_executable(configcp):
  try:
    print '--- Check that the executable ('+ configcp.get("main","executable")  +')is present in '+path
    f = open(configcp("main", "executable"), 'r')
    f.close()
  except:
    print '### Can not find ' + configcp("main", "executable")
    sys.exit()
  print '... executable found. Going ahead'
	



def parse_arguments():
  print '--- Parsing user arguments'
  parser = OptionParser()
  parser.add_option( "--config-file")
  
  parser.add_option("--search",
      dest='search',default='BNS',
      help=" <BNS, BBH, PBH , BHNS, S5 (1,60)>")
  parser.add_option("--bank-ffinal",
      dest='bank_ffinal', default=2047, type='float',
      help="upper frequency to be used" )
  parser.add_option("--max-total-mass",
      dest='max_total_mass', default=-1, type='float',
      help="max total mass (injection)" )
  parser.add_option("--fast-simulation",
      action="store_true", default="false",
      dest='fast_simulation', 
      help="fast simulation option" )
  parser.add_option("--bhns-injection",
      action="store_true", default="false",
      dest='bhns_injection', 
      help="bhns injection only. If search arguments is set to BHNS, this parameter will always be used." )



  (options, args) = parser.parse_args()
  return options,args

# -----------------------------------------------------------------------------------
options, args = parse_arguments()
configcp = ConfigParser.ConfigParser()
configcp.read(options.config_file)
        
    
arguments = create_condor_file(configcp)
print """
	The condor script will use the following arguments 
	-------------------------------------------
    """
print arguments
print '\n--- The number of simulation requested is '+configcp.get("simulation", "ntrial")
print '--- They will be split into '+ configcp.get("simulation", "njobs")+' jobs'

# create the condor file using the input parameter stored in BE
create_finalise_script(configcp)
create_dag_file(configcp)
create_bank(configcp, arguments)

print '--- Generating the prototype xml file for merging condor job'
arguments = ""
command = configcp.get("main", "executable") + ' ' + arguments +' --print-prototype 1>bep_proto.out 2>bep_proto.err'
os.system(command)
print '... done'
time.sleep(.5)
    
print """--- In order to start the job, type
--------------------------------------------
condor_submit_dag -maxjobs 100  bep.dag

or 

condor_submit_dag -maxjobs 100  -f bep.dag
--------------------------------------------
                
Once the dag is finished and all the job are completed, get back all
the results together within an xml file by using the script called : finalise.sh

Ideally, this script should be put within the daga file"""

create_finalise_script(configcp)
os.system('mkdir log')
        
#if __name__ == "__main__":
#    main()
