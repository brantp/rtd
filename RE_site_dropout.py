#!/usr/bin/env python

from scipy.stats.distributions import norm,poisson
import sys,numpy,time

def count_likelihood_standard(this_counts,tot_counts,this_port,this_cv,this_fract_recov=1):
    L = [norm.pdf(this_c,loc=tot_c*this_port*this_fract_recov,scale=(tot_c*this_port*this_fract_recov)*this_cv) \
         for this_c,tot_c in zip(this_counts,tot_counts)]
    logL = numpy.log2(numpy.array(L))
    return sum(logL)

# <codecell>

def count_likelihood_poisson(this_counts,tot_counts,this_port,this_cv,this_fract_recov=1):
    L = [poisson.pmf(this_c,tot_c*this_port*this_fract_recov) \
         for this_c,tot_c in zip(this_counts,tot_counts)]
    logL = numpy.log2(numpy.array(L))
    return sum(logL)

# <codecell>

def fake_null_cuts(this_counts,freq_cut_null):
    num_absent = int(round(freq_cut_null**2 * len(this_counts)))
    num_het = int(round((2 * freq_cut_null * (1 - freq_cut_null)) * len(this_counts)))
    ret = []
    for i in range(num_absent):
        ret.append(0)
    for i in range(num_absent,num_absent+num_het):
        ret.append(int(this_counts[i]*numpy.random.normal(0.5,0.01)))
    for i in range(num_absent+num_het,len(this_counts)):
        ret.append(this_counts[i])
    return numpy.array(ret)

# <codecell>

def count_likelihood_skew(this_counts,tot_counts,this_port,this_cv,freq_cut_null,err=1e-5,Lfn=count_likelihood_standard):
    num_absent = int(round(freq_cut_null**2 * len(this_counts)))
    num_het = int(round((2 * freq_cut_null * (1 - freq_cut_null)) * len(this_counts)))
    ind_absent = [i for i,ncounts in sorted(enumerate(this_counts),key = lambda x:x[1])[:num_absent]]
    ind_het = []
    for i, portion in sorted(enumerate(this_counts/tot_counts.astype('float')),key = lambda x:x[1]):
        if len(ind_het) == num_het:
            break
        if i in ind_absent: 
            continue
        ind_het.append(i)
    ind_full = list((set(range(len(this_counts)))-set(ind_het))-set(ind_absent))
    #print [(i,this_counts[i],(this_counts/tot_counts.astype('float'))[i]) for i in ind_absent]
    #print [(i,this_counts[i],(this_counts/tot_counts.astype('float'))[i]) for i in ind_het]
    #print [(i,this_counts[i],(this_counts/tot_counts.astype('float'))[i]) for i in ind_full]
    #print 'absent %s\nhet %s\nfull %s' % (ind_absent,ind_het,ind_full)
    L0li = [err**this_counts[i] for i in ind_absent]
    L1 = Lfn(this_counts[ind_het],tot_counts[ind_het],this_port,this_cv,0.5)
    L2 = Lfn(this_counts[ind_full],tot_counts[ind_full],this_port,this_cv,1)
    return sum(numpy.log2(L0li)) + L1 + L2

# <codecell>

def max_likelihood_standard(this_counts,tot_counts, \
                            num_portions_to_test=10,mincv=0.1,maxcv=2.0,num_cv_to_test=10, \
                            npass=2,take_top=0.1,Lfn=count_likelihood_standard):
    this_portions = this_counts/tot_counts.astype('float')
    
    test_ports = numpy.linspace(this_portions.min(),this_portions.max(),num_portions_to_test)
    if Lfn == count_likelihood_poisson:
        test_cvs = [1]
    else:
        test_cvs = numpy.linspace(mincv,maxcv,num_cv_to_test)
    #print this_portions,test_ports,test_cvs
    for i in range(npass):
        l_out = sorted([ (Lfn(this_counts,tot_counts,test_port,test_cv),\
                  (test_port,test_cv)) \
                for test_port in test_ports for test_cv in test_cvs])
        top = int(len(test_ports)*len(test_cvs)*take_top)
        #print >> sys.stderr, l_out[-1]
        top_ports,top_cv = zip(*[el[1] for el in l_out[-top:]])
        test_ports = numpy.linspace(min(top_ports),max(top_ports),num_portions_to_test)
        test_cvs = numpy.linspace(min(top_cv),max(top_cv),num_cv_to_test)
    return l_out[-1]

# <codecell>

def max_likelihood_skew(this_counts,tot_counts, \
                        num_portions_to_test=10,mincv=0.1,maxcv=2.0,num_cv_to_test=10, \
                        num_freq_to_test=10,npass=2,take_top=0.1,Lfn=count_likelihood_standard):
    this_portions = this_counts/tot_counts.astype('float')
    test_ports = numpy.linspace(this_portions.min(),this_portions.max(),num_portions_to_test)
    if Lfn == count_likelihood_poisson:
        test_cvs = [1]
    else:
        test_cvs = numpy.linspace(mincv,maxcv,num_cv_to_test)
    test_NFs = numpy.linspace(0,1,num_freq_to_test)
    for i in range(npass):
        l_out = sorted([ (count_likelihood_skew(this_counts,tot_counts,test_port,test_cv,test_NF,Lfn=Lfn),\
                  (test_port,test_cv,test_NF)) \
                for test_port in test_ports for test_cv in test_cvs for test_NF in test_NFs])
        top = int(len(test_ports)*len(test_cvs)*take_top)
        #print >> sys.stderr, l_out[-1]
        top_ports,top_cv,top_NF = zip(*[el[1] for el in l_out[-top:]])
        test_ports = numpy.linspace(min(top_ports),max(top_ports),num_portions_to_test)
        test_cvs = numpy.linspace(min(top_cv),max(top_cv),num_cv_to_test)
        test_NFs = numpy.linspace(min(top_NF),max(top_NF),num_cv_to_test)
    return l_out[-1]   
    

def run_main(tot_mean,tot_sd,n_ind,this_mean,this_cv,this_NF,n_steps,n_pass,Lfn=count_likelihood_standard):
    
    #tot_mean = 300000
    #tot_sd = 150000
    #n_ind = 48
    tot_counts = numpy.array([x>10000 and x or 10000 \
                              for x in numpy.random.normal(tot_mean,tot_sd,size=n_ind).astype('int')])

    #this_mean = 10
    port = this_mean/float(tot_mean)
    #this_cv = 0.5
    this_counts = numpy.array([x>0 and x or 0 for x in tot_counts*numpy.random.normal(port,port*this_cv,n_ind)]).astype('int')

    print >> sys.stderr, tot_counts,numpy.mean(tot_counts),numpy.std(tot_counts)
    print >> sys.stderr, this_counts,numpy.mean(this_counts),numpy.std(this_counts)
    mean_port =  numpy.mean(this_counts/tot_counts.astype('float'))
    std_port = numpy.std(this_counts/tot_counts.astype('float'))
    print >> sys.stderr, mean_port, std_port/mean_port

    L0,(L0_portion,L0_cv) = max_likelihood_standard(this_counts,tot_counts,num_portions_to_test=n_steps,num_cv_to_test=n_steps,npass=n_pass,Lfn=Lfn)
    L1,(L1_portion,L1_cv,F) = max_likelihood_skew(this_counts,tot_counts,num_portions_to_test=n_steps,num_cv_to_test=n_steps,num_freq_to_test=n_steps,npass=n_pass,Lfn=Lfn)
    print >> sys.stderr, 2*(-L0 - -L1), L0_portion,L0_cv
    LRTneg = 2*(-L0 - -L1)

    skew_counts = fake_null_cuts(this_counts,this_NF)
    L0,(L0_portion,L0_cv) =  max_likelihood_standard(skew_counts,tot_counts,num_portions_to_test=n_steps,num_cv_to_test=n_steps,npass=n_pass,Lfn=Lfn)
    L1,(L1_portion,L1_cv,F) =  max_likelihood_skew(skew_counts,tot_counts,num_portions_to_test=n_steps,num_cv_to_test=n_steps,num_freq_to_test=n_steps,npass=n_pass,Lfn=Lfn)
    print >> sys.stderr, 2*(-L0 - -L1), L1_portion,L1_cv,F
    LRTpos = 2*(-L0 - -L1)

    return {'LRT':{'neg':LRTneg,'pos':LRTpos},'F':F}

if __name__ == "__main__":

    tot_mean,tot_sd,n_ind,this_mean,this_cv,this_NF,n_steps,n_pass = map(float,sys.argv[1:])
    n_ind, n_steps, n_pass = map(int,[n_ind, n_steps, n_pass])

    print run_main(tot_mean,tot_sd,n_ind,this_mean,this_cv,this_NF,n_steps,n_pass)
