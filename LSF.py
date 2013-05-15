'''
#####
# generic LSF functions
#####

core functionality resides in lsf_jobs_submit and lsf_wait_for_jobs
see esp. lsf_jobs_submit doctring for usage
'''


import subprocess,os,sys,re,time,numpy
import random, string

USER = os.environ['USER']


def random_filename(chars=string.hexdigits, length=8, prefix='', suffix='', \
                        verify=True, attempts=10):
    for attempt in range(attempts):
        filename = ''.join([random.choice(chars) for i in range(length)])
        filename = prefix + filename + suffix
        if not verify or not os.path.exists(filename):
            return filename

def clear_lsf_scratch(username=USER):
    pass

def lsf_jobs_dict(username=USER):
    '''
    returns status dictionary of currently running jobs
    currently horks if parallel jobs are running (since there are lines with only one item)
    fix w/ regex - done (not yet tested)

    better yet, replace with a wrapper around lsf_get_job_details.
    (or replace all uses with lsf_get_job_details calls)
    '''
    
    statpairs = []
    for l in subprocess.Popen(['bjobs','-w'],stdout=subprocess.PIPE,stderr=open('/dev/null','w')).stdout.readlines()[1:]:
        match = re.search('^(\d+)\s+%s\s+(.+?)\s' % username,l)
        if match:
            statpairs.append(match.groups())
    return dict(statpairs)

def lsf_get_run_time(jids):
	'''given a list of job ids, returns two dicts
	the first contains jid:time pairs for completed jobs, the second for incomplete jobs
	
	given a single job id, returns the time that job ran or has run for'''

	import time

	def jobtime_to_secs(timestr):
		return time.mktime(time.localtime()[0:1]+tuple([int(s) for s in re.split('[:/-]',timestr)])+time.localtime()[-3:])
	
	def get_time(jid):
		try:
			start,end = subprocess.Popen('bjobs -W %s 2> /dev/null' % jid,shell=True,stdout=subprocess.PIPE).stdout.readlines()[1].split()[-2:]
		except IndexError:
			return ('GONE',None)
		if end == '-':
			if start == '-':
				return ('PEND',None)
			else:
				start = jobtime_to_secs(start)
				end = time.mktime(time.localtime())
				return ('CURR',end-start)
		else:
			start = jobtime_to_secs(start)
			end = jobtime_to_secs(end)
			return ('DONE',end-start)
	
	if isinstance(jids,str):
		return get_time(jids)
	else:
		times = {}
		times['PEND']={}
		times['CURR']={}
		times['DONE']={}
		times['GONE']={}
		for j in jids:
			stat,secs = get_time(j)
			times[stat][j] = secs
		return times['DONE'],times['CURR']
	

def lsf_get_job_details(jid):
	stat = subprocess.Popen('bjobs -l %s 2> /dev/null' % jid, shell=True, stdout=subprocess.PIPE).stdout.read()
	statlines = []
	for i in stat.replace(' >\n','>\n').replace('>\n','>, ').split('\n'):
		if i.startswith(' '*21):
			i = i[21:]
		statlines.append(i.rstrip('\n'))
	statstr = ''.join(statlines)
    #statstr = ''.join([i.lstrip().rstrip('\n') for i in stat.replace('>\n','>, ').split('\n')])
    #statstr = statstr.split(';')[0]
	statstr = re.split(r';[^<>]+Started\son',statstr)[0]
        #print >> sys.stderr, statstr
	statli = re.findall('\s*?(.+?)\s?<(.+?[^\s])>,?',statstr.replace(' >,','>,'))
	statdict = dict([(i.strip(),j.strip()) for i,j in statli])
        if statdict.get('Command',None):
		statdict['Command'] = statdict['Command'].replace('"',r'\"') #must re-escape command string
	return statdict

def lsf_jobs_submit(cmds,outfile,queue='short_serial',bsub_flags='',jobname_base=None, prereq_names=None, num_batches=None, batch_size=None, sh_tmp=os.path.expanduser('~')+'/tmp/sh'):
    '''given a list of bash-runnable commands (strings), submits to indicated queue, tags each submission with <bsub_flags>
	(-R resource requirements, e.g.)

	if jobname_base supplied, returns (JOBIDS,NAMEDICT) tuple, and names each job <jobname_base><random_id>
	
	batching is optional; if num_batches is supplied, overrides batch_size with len(cmds)/num_batches
	
	if batches are used, or if commands exceed 1024 characters, .sh scripts will be generated, placed in sh_tmp and submitted to bsub in lieu of strings in cmds

	logs will be written to outfile (i.e. -o <outfile>)

	usage ex:
	#some code to generate a list of commands (e.g. that you could run in a bash shell)
	cmds = ["ls | wc > /tmp/somefile"] * 20
	#then run them out, in this case 20 commands will end up in 10 batches, 2 jobs per, so make sure <sh_tmp> is writable for .sh script creation
	jobids,namedict = LSF.lsf_jobs_submit(cmds,outfile="/tmp/logfile",queue="normal_serial",jobname_base="ls-wordcount",num_batches=10,bsub_flags="-R \"select[mem > 15000]\"")
	
	#see docstr for lsf_wait_for_jobs for details.  
	#supplying restart_outfile invokes SSUSP checking.
	#supplying restart_z invokes slow job scrubbing (i.e. for jobs that land on overworked nodes).
	#DONT USE restart_z WITH HETEROGENEOUS RUN DURATION JOBS!
	LSF.lsf_wait_for_jobs(jobids,restart_outfile="/tmp/restarts-log",namedict=namedict,restart_z=6)
    '''
    if isinstance(cmds,str):
        cmds = [cmds]

    #print >> sys.stderr, '\n\n',num_batches,'\n\n'

    if num_batches:
        tot = len(cmds)
        batch_size = tot / num_batches
        print >> sys.stderr, 'batching invoked, %s batches requested (%s jobs per batch)' % (num_batches,batch_size)

    if batch_size:
        cmds = ['; '.join(cmds[i:i+batch_size]) for i in xrange(0,len(cmds),batch_size)]

    try:
        os.makedirs(os.path.dirname(outfile))
    except:
        pass
    
    jobids = {}
    if jobname_base:
        namedict = {}
    else:
        namedict = None
    execbase = 'bsub -q %s %s' % (queue,bsub_flags)
    print >> sys.stderr,'Adding jobs'

    for execstr in cmds:
        print >> sys.stderr,'.',

        time.sleep(0.1) # try not to stompy LSF

        execfull = execbase
        if jobname_base:
            jname = random_filename(prefix=jobname_base)
            execfull += ' -J %s -o %s-%s.lsflog' % (jname,outfile,jname)
        else:
            jname = None
            execfull += ' -o %s.lsflog' % random_filename(prefix='noName-')
            
        if ';' in execstr or len(execstr) > 1024 or batch_size:
            if jname:
                sh_tmp_path = os.path.join(sh_tmp,jobname_base)
                sh_name = os.path.join(sh_tmp_path,jname+'.sh')
            else:
                sh_tmp_path = sh_tmp
                sh_name = os.path.join(sh_tmp_path,random_filename()+'.sh')
            try:
                os.makedirs(sh_tmp_path)
            except OSError:
                pass

            open(sh_name,'w').write('#!/usr/bin/env sh\nset -e\n'+execstr.replace('\\"','"'))
            os.system('chmod +x '+sh_name)
            execstr = sh_name
            
        if prereq_names:
            if isinstance(prereq_names,str):
                prereq_names = [prereq_names]
            execfull += ' -w "%s"' % '&&'.join(['done(%s)' % pn for pn in prereq_names])

        match = None
        while match is None:
           ret = subprocess.Popen(execfull + ' "%s"' % (execstr),shell=True,stdout=subprocess.PIPE).stdout.read()
           match = re.search('<(\d+)>',ret) 

        jid = match.groups()[0]
        jobids[jid] = execstr
        if jobname_base:
            namedict[jid] = jname
            
    if namedict is not None:
        return jobids,namedict
    else:
        return jobids

def lsf_jobs_status(jobids):
	if not jobids:
		return {}
	alljobs = lsf_jobs_dict()
	jstats = {}.fromkeys(set([v for k,v in alljobs.items() if k in jobids.keys()]+['DONE']),0)
	for j in jobids.keys():
		try:
			jstats[alljobs[j]] += 1
		except KeyError:
			jstats['DONE'] += 1
	return jstats

def lsf_wait_for_jobs(jobids,restart_outfile=None,restart_queue='normal_serial',sleeptime = 20,namedict=None,restart_z=None,restart_stragglers_after=0.75,kill_if_all_ssusp=False):
	'''waits until all listed jobs are finished.  
	specifying a restart_outfile (log for restarted items) will invoke checking for and re-starting SSUSP jobs
	(which is HIGHLY recommended)
	
	moderately experimental:
	if restart_z is not None (i.e. is a number) and fraction completed > <restart_stragglers_after>
	check for jobs w/ runtime > <restart_z> stdevs above mean run for completed jobs
	NB 20090908 this behavior is hardcoded to avoid invocation until jobtime exceeds 30sec (line 228)
	
	embarassing: 
	pretty sure restart_queue doesnt do anything anymore.  sad.
	'''
	
	import time,datetime
	print >> sys.stderr, 'running %s jobs' % len(jobids)
	status = lsf_jobs_status(jobids)
        t = time.time()
        maxllen = 0
	while any([status[k] for k in status.keys() if k != 'DONE']):
		time.sleep(sleeptime)
		status = lsf_jobs_status(jobids)
		pctdone = status['DONE'] / float(sum(status.values()))
                if kill_if_all_ssusp and set(status.keys()) == set(['DONE','SSUSP']):
			sjobs = [j for j,s in lsf_jobs_dict().items() if s == 'SSUSP' and j in jobids.keys()]
                        print >> sys.stderr, 'all remaining jobs %s in SSUSP; terminating' % sjobs
                        lsf_kill_jobs(sjobs)
                        continue
		if 'SSUSP' in status.keys() and restart_outfile:
			print >> sys.stderr,'restarting',status['SSUSP']
                        if namedict:
                            jobids,namedict = lsf_restart_susp(jobids,restart_outfile,restart_queue,namedict=namedict,remove_old=True)
                        else:
                            jobids = lsf_restart_susp(jobids,restart_outfile,restart_queue,remove_old=True)
		if restart_z and pctdone > restart_stragglers_after:
			done,curr = lsf_get_run_time(jobids.keys())
			runtimes = numpy.array(done.values())
                        outl = '\r%s %s (%3d%% av %0.1f st %0.1f sec)' % (str(datetime.timedelta(seconds=int(time.time() - t))),status.__repr__(),pctdone*100,runtimes.mean(),runtimes.std())
                        if len(outl) < maxllen:
                            pad = maxllen - len(outl)
                            outl += ' '*pad
                        else:
                            maxllen = len(outl)
			sys.stderr.write(outl)
                        sys.stderr.flush()
                        
			cut = runtimes.mean() + (restart_z * runtimes.std())
			restarts = []
			for j,t in curr.items():
				if t > cut and t > 30: #hardcoded 30sec floor on kills
					print >> sys.stderr, 'job %s has run for %0.1f sec, cutoff of %0.1f set based on run mean/stdev of %0.1f/%0.1f' % (j,t,cut,runtimes.mean(),runtimes.std())
					restarts.append(j)
			if restarts:
				restartids, restartnames, kills = lsf_restart_jobs(restarts,return_kills=True)
				for j in kills:
					del jobids[j]
				jobids.update(restartids)
				namedict.update(restartnames)
		else:
			outl = '\r%s %s (%3d%%)' % (str(datetime.timedelta(seconds=int(time.time() - t))),status.__repr__(),pctdone*100)
                        if len(outl) < maxllen:
                            pad = maxllen - len(outl)
                            outl += ' '*pad
                        else:
                            maxllen = len(outl)
                        
			sys.stderr.write(outl)
                        sys.stderr.flush()
        print >> sys.stderr, '\ncompleted iteration in',str(datetime.timedelta(seconds=int(time.time() - t)))

def lsf_run_until_done(to_run_dict,logfile,queue,bsub_flags,jobname_base,num_batches,MAX_RETRY):
    from run_safe import unfinished_cmds
    cmds = unfinished_cmds(to_run_dict)

    retries = 0
    last_cmds = []
    while len(cmds) > 0:
        print >> sys.stderr, '%s: %s cmds to run in %s batches on queue %s, logs in %s' % (jobname_base,len(cmds),num_batches,queue,logfile)
        #code to halt execution on recurrent errors
        if set(last_cmds) == set(cmds):
            if retries > MAX_RETRY:
                errstr = 'maximum number of retry attempts (%s) exceeded with identical jobs lists.  Check logs (%s) for recurrent errors' % (MAX_RETRY,logfile)
                raise IOError, errstr
            else:
                retries += 1
        last_cmds = cmds

        jobids,namedict = lsf_jobs_submit(cmds,logfile,queue,bsub_flags,jobname_base=jobname_base,num_batches=num_batches)
        time.sleep(20)
        lsf_wait_for_jobs(jobids,logfile,namedict=namedict)

        cmds = unfinished_cmds(to_run_dict)
    print >> sys.stderr, 'DONE\n'


def lsf_kill_jobs(jobids):
    if isinstance(jobids,dict):
        kills = jobids.keys()
    else:
        kills = jobids
    if not isinstance(kills,list):
        kills = [kills]
    execstr = 'bkill -r  %s > /dev/null 2> /dev/null' % ' '.join(kills)
    #print >> sys.stderr, execstr
    os.system(execstr)
    return kills

def lsf_restart_jobs(jobids,keep_prereqs=True,return_kills=False):
    if isinstance(jobids,dict):
        kills = jobids.keys()
    else:
        kills = jobids
    if not isinstance(kills,list):
        kills = [kills]

    restartids = {}
    restartnames = {}    

    kills_finish = []
    kills_skip = []
    for j in kills:
        jdet = lsf_get_job_details(j)
        depstr = jdet.get('Dependency Condition',None)
        if depstr and keep_prereqs:
            deps = re.findall('done\((.+?)\)',depstr)
        else:
            deps = None
        try:
            jids,ndict = lsf_jobs_submit([jdet['Command']],jdet.get('Output File','/tmp/%s_restarted_jobs' % USER),queue=jdet['Queue'],jobname_base=jdet['Job Name']+'-', prereq_names=deps)
            kills_finish.append(j)
            restartids.update(jids)
            restartnames.update(ndict)
        except KeyError:
            kills_skip.append(j)
    if kills_skip:
        print >> sys.stderr, 'no restart data available for %s; not killed' % kills_skip
    lsf_kill_jobs(kills_finish)

    if return_kills:
        return restartids, restartnames, kills_finish
    else:
        return restartids, restartnames
	
def lsf_restart_susp(jobids,outfile='/tmp/restarted_jobs',queue='normal_serial',return_only_restarted_ids=False,namedict=None,remove_old=False):
	'''current logic just restarts the job as-is, from lsf_get_job_details()
	obsolete: N.B. if namedict supplied, will return (JOBID, NAMEDICT) tuple!'''
	stat = lsf_jobs_dict()
	susp = [k for k,v in stat.items() if v == 'SSUSP' and k in jobids.keys()]
	if any(susp):
		restartids, restartnames, kills = lsf_restart_jobs(susp,return_kills=True)
		if return_only_restarted_ids:
			if namedict:
				return restartids,restartnames
			else:
				return restartids
		else:
			if remove_old:
				for k in kills:
					try:
						del jobids[k]
						#for now, retain killed jobs in namedict; this could change	
					except KeyError:
						pass

			if namedict:
				jobids.update(restartids)
				namedict.update(restartnames)
				return jobids,namedict
			else:
				jobids.update(restartids)
				return jobids

def lsf_no_success_from_log(logfile,success_funct=any):
	'''returns a list of commands that never finished successfully in the supplied lsf log
	
	uses success_funct to determine "passing" command strings.
	consider:
	success_func=lambda x:x[-1]
	to check that the _last_ time a command was attempted, it completed successfully'''
	
	log = open(logfile).read()
	results = re.findall('# LSBATCH: User input\\n(.+?)\n([^#]+?)output',log)
	d = {}.fromkeys(set([i[0] for i in results]),None)
	for k in d:
		d[k] = []
	for cmd,result in results:
		d[cmd].append('Success' in result)
		
	return [k for k,v in d.items() if not success_funct(v)]
	
def lsf_success_from_log(logfile,success_funct=any):
	'''essentially returns the inverse of lsf_no_success_from_log (big surprise)'''
	
	log = open(logfile).read()
	results = re.findall('# LSBATCH: User input\\n(.+?)\n([^#]+?)output',log)
	d = {}.fromkeys(set([i[0] for i in results]),None)
	for k in d:
		d[k] = []
	for cmd,result in results:
		d[cmd].append('Success' in result)
		
	return [k for k,v in d.items() if success_funct(v)]
	
