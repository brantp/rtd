#!/usr/bin/env python

import os,sys

def safe_script(cmd,donefile,donesuffix='.done',scriptfile=None,force_write=False):
    if scriptfile is None:
        scriptfile = donefile+'.sh'
    if os.path.exists(scriptfile) and not force_write:
        #file present; use
        pass
    else:
        try:
            os.unlink(scriptfile)
        except:
            pass
        open(scriptfile,'w').write('#!/usr/bin/env sh\nset -e\n%s\n' % (cmd))
        ret = os.system('chmod +x %s' % scriptfile)
        if ret != 0:
            raise OSError, 'cannot set execute for %s' % scriptfile
    return 'run_safe.py \"%s\" %s' % (scriptfile,donefile+donesuffix)
                               
def unfinished_cmds(to_run_dict,finished_ext='.done'):
    cmds = []
    for finished_base,cmd in to_run_dict.items():
        if not os.path.exists(finished_base+finished_ext):
            cmds.append(cmd)
    return cmds

if __name__ == "__main__":
    cmd,done = sys.argv[1:]

    if os.path.exists(done):
        print >> sys.stderr, 'completion flag exists: %s' % done
    else:
        ret = os.system(cmd)
        if ret == 0:
            os.system('touch %s' % done)

            
