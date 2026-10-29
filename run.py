
import numpy as np
import matplotlib.pyplot as plt
import os
import subprocess
import time
import scipy.stats as st
from scipy.optimize import curve_fit
from scipy.special import erf
from scipy.stats import sem
import argparse

#This is the original code used in the thesis, I might update this later to be more clear, but I recommend making your own to use with the program
#----------------------------------------------------------------
#   Input parameters
#
#
#----------------------------------------------------------------


parser= argparse.ArgumentParser(add_help=False)
parser.add_argument("-t","--transition",help="The type of transition to use, current options: fast,slow,equal_size",required=True)
parser.add_argument("--input",help="folder to data to plot instead of running simulation",type=str,default="")
parser.add_argument("--plot",help="Whether to plot or not",type=bool)
args,_ = parser.parse_known_args()
read_data=args.input
do_plot=args.plot
print(do_plot)
ttype=args.transition.lower()
if not read_data:
    if ttype=="fast":
        parser.add_argument("-G_f",help="value of percolation rate Gamma(t_f)=Gamma_f",required=True)
        #parser.add_argument("-tau_f",help="maximum value of tau parameter, tau=beta*(t-t_f)",required=True)
    elif ttype=="slow":
        parser.add_argument("-G_0",help="value of constant percolation rate Gamma_0 for vacuum",required=True)
        parser.add_argument("-eta_f",help="maximum value of conformal time eta",required=True)
        parser.add_argument("-h0","--hubble",help="value of the hubble constant at start of transition",default=1)
        parser.add_argument("-N_f","--N_final",help="",required=True);
    elif ttype=="equal_size":
        parser.add_argument("-r","--radius",help="bubble radii",required=True)
    else:
        raise ValueError(f"Unknown transition type: {args.transition}\n Available types are (fast,slow,equal_size)")

    parser.add_argument("-j","--threads",help="number of threads to use, total trials is threads*n_repeat",required=True,type=int)
    parser.add_argument("-n","--n_repeat",help="number of trials each thread does",required=True,type=int)
    parser.add_argument("-b","--boundaries", help="list of boundary sizes to simulate",nargs='+',required=True,type=int)
    parser.add_argument("--step",help="list of boundary sizes to simulate with step",nargs=1,required=False,type=int)
    

parser.add_argument(
"-h", "--help",
action="help",
default=argparse.SUPPRESS,
help="show this help message and exit"
)
args=parser.parse_args()


#Main repeats

if not read_data:
    repeats=args.n_repeat
    threads=args.threads
    step=args.step

    simulation_name=""
    if ttype=="slow":
        p_f=args.G_0
        simulation_name="slow_transition_sim.o"
    elif ttype=="fast":
        p_f=args.G_f
        simulation_name="fast_transition_sim.o"
    elif ttype=="equal_size":
        p_f=args.radius
        simulation_name="equal_sized_spheres.o"    

    boundaries=args.boundaries
    if step:
        boundaries=np.arange(boundaries[0],boundaries[1],step[0])
#----------------------------------------------------------------
#   Definitions
#
#
#----------------------------------------------------------------

#def step_function(tau,tau_critical,nabla):
#    return 1/2*(1+np.tanh((tau-tau_critical)/nabla))
def step_function(tau,tau_critical,nabla):
    return 1/2*(1+erf((tau-tau_critical)/nabla))
def wilson(c, p_hat,n):
    z_score = st.norm.ppf(1 - ((1 - c) / 2))
    high=(1+z_score**2/n)**(-1)*(p_hat+z_score**2/(2*n)+z_score/(2*n)*np.sqrt(4*n*p_hat*(1-p_hat)+z_score**2))
    low=   (1+z_score**2/n)**(-1)*(p_hat+z_score**2/(2*n)-z_score/(2*n)*np.sqrt(4*n*p_hat*(1-p_hat)+z_score**2)) 
    return p_hat-low,high-p_hat
def wilson2(c, p_hat,n):
    z_score = st.norm.ppf(1 - ((1 - c) / 2))
    high=(1+z_score**2/n)**(-1)*(p_hat+z_score**2/(2*n)+z_score/(2*n)*np.sqrt(4*n*p_hat*(1-p_hat)+z_score**2))
    low=   (1+z_score**2/n)**(-1)*(p_hat+z_score**2/(2*n)-z_score/(2*n)*np.sqrt(4*n*p_hat*(1-p_hat)+z_score**2)) 
    return low,high
def wald(c,p_hat,n):
	z_score = st.norm.ppf(1-((1-c)/2))
	z_score/np.sqrt(n)*np.sqrt(p_hat*(1-p_hat))
	return  z_score/np.sqrt(n)*np.sqrt(p_hat*(1-p_hat))
def finitescaling(L,phi_c,nu,a):
    return L**(-1/nu)*(a)+phi_c

#----------------------------------------------------------------
#   File manipulation, the data gets saved to /tmp/data.dat momentarily
#   if this folder already exists throws exception
#----------------------------------------------------------------
if not read_data:
    try:
        os.mkdir("tmp")
    except FileExistsError:
        answ=input("tmp exists, remove? (y/n)")
        
        if answ.lower()=="y":
            os.system("rm -r ./tmp/")
            os.mkdir("tmp")
            print("Removed tmp, running...")
        else:
            print("Quitting")
            exit()
#----------------------------------------------------------------
#   Main loop
#
#
#----------------------------------------------------------------

#For initial value forking
    repeat_step=2
    repeat_fork=4



    if ttype=="equal_size":
        tau_div=14
        tau_upper=0.290
        tau_lower=0.287
        adjust_size=0.001
    else:
        tau_div=10
#Initial starting points to refine from
        tau_upper=-1
        tau_lower=-1.04
        adjust_size=0.05
#Aim for 0.4-0.6 for symmetric error at 1000 repeats, i.e must have 0.4 or 0.5

#Do quick test to refine this



    start=time.time()
    print(tau_lower,tau_upper)
#NOTE: unnecessary for equal_sized_spheres, need to adjust for slow/fast ones
    for _ in range(repeat_fork):
        taus = np.linspace(tau_lower,tau_upper,tau_div)
        for tau in taus:
            os.system(f"./{simulation_name} {tau} {p_f} {boundaries[0]} {repeat_step} {threads} >> ./tmp/data_test.dat")
            # print(os.system("cat ./tmp/data_test.dat"))
        taus, successes = np.loadtxt(f"./tmp/data_test.dat",unpack=True)
        for i,(tau,success) in enumerate(zip(taus,successes)):
            if success>successes[i-1]:
                prob=success/(repeat_step*threads)
                # print("debug",prob,tau,tau_lower,tau_upper)
                if ( prob<0.3 and tau>tau_lower):
                    tau_lower=tau
                elif ( prob>0.6 and tau<tau_upper):
                    tau_upper=tau

        #probabilities = successes/(threads*repeats)
        print(tau_lower,tau_upper)

        os.system("rm ./tmp/data_test.dat")

    print(f"Refinement {time.time()-start}s")



    taus = np.linspace((tau_lower+tau_upper)/2,1,tau_div)
    tau=tau_lower-(tau_upper-tau_lower)/3
    c=["r","g","b","orange","black"]
    i=-1
    tau_crit=[]
    tau_crit_err=[]
    plt.figure(1)
    for boundary in boundaries:
        i+=1
        os.system("touch ./tmp/data.dat")
        print("boundary,tau:",boundary,tau)
        start=time.time()
        for _ in range(tau_div):
            os.system(f"./{simulation_name} {tau} {p_f} {boundary} {repeats} {threads} >> ./tmp/data.dat")
            taus, successes = np.loadtxt(f"./tmp/data.dat",unpack=True)
            try:
                probability = successes[-1]/(repeats*threads) 
            except IndexError:
                probability=successes/(repeats*threads)
            if (probability<0.35): 
                tau=tau+adjust_size #0.03 initially
            elif (probability >= 0.65):
                tau=tau-adjust_size*np.random.random()
            elif (probability >=0.35 ):
                tau=tau+adjust_size/3

        print(f"{(time.time()-start)/60}min")

        #----------------------------------------------------------------
        #   Data manipulation
        #
        #
        #----------------------------------------------------------------



        error_both=[]
        error_max=[]
        taus,successes = np.loadtxt(f"./tmp/data.dat",unpack=True)
        probabilities = successes/(threads*repeats)
        for tau_l,probability in zip(taus,probabilities):
            error = wilson(0.69,probability,threads*repeats)
#        error_max.append(abs(error[1]-error[0])/(2*1.96))
            error_max.append(max(error))
            error_both.append(error)
            if probability>=0.35 and tau_l<tau:
                tau=tau_l
            print(tau_l,probability,error)
        data=np.array([taus,probabilities])
        data=data[:,np.argsort(data[0])]
        mask=np.logical_and(data[1]>=0.35,data[1]<=0.65) #Mask out data not in correct range
        taus=data[0][mask]
        probabilities=data[1][mask]
        error_max=np.array(error_max)[mask]
        error_both=np.array(error_both)[mask]
        print(probabilities)
        print(successes)
        #TODO: Handle failure so we dont crash here 
        par,cov=curve_fit(step_function,taus,probabilities,p0=[-0.2,1],sigma=error_max,absolute_sigma=True)
        fit_error = np.sqrt(np.diag(cov))
        error_both=np.array(error_both)

        tau_crit_err.append(fit_error[0])
        tau_crit.append(par[0])

        data=np.array(np.hstack([par,fit_error]))
        data = np.vstack([data,np.vstack([taus,probabilities,error_both[:,0],error_both[:,1]]).T])

        plt.errorbar(taus,probabilities,yerr=error_max,linestyle="",capsize=4,c=c[i%len(c)],label=f"L={boundary}")
        plt.scatter(taus,probabilities,c=c[i%len(c)])

        x=np.linspace(0,1,1000)
        plt.plot(x,step_function(x,*par),c=c[i%len(c)],label=f"fit L={boundary}")
        
        comment = ("Data is stored in the following manner\n\
        First row:         tau_crit,nabla,tau_crit_fit_error,nabla_fit_error\n\
        Rows after this:   tau,probability, wilson_low, wilson_high")
        current_time=time.strftime('%d-%h--%H:%M',time.localtime(time.time()))
        np.savetxt(f"data_{boundary}_{current_time}_wil_{simulation_name}.dat",data,header=f"Boundary size = {boundary},repeat={repeats},threads={threads}\n {comment}")
        os.system("rm ./tmp/data.dat")
    ax =plt.gca()
    ax.set_xlabel(r"\tau")
    ax.set_ylabel(r"$Probability$")
    plt.legend()
    plt.savefig(fname=f"fig_1_{current_time}_{ttype}.png")
    plt.figure(2)

if read_data:
    
    def system_call(cmd,live_output=False):
        with subprocess.Popen(cmd.split(" "),stdout=subprocess.PIPE,stderr=subprocess.PIPE) as proc:
            if live_output:
                for line in proc.stdout:
                    print(str(line,'utf-8').rstrip('\n')) #Note that for example pip's progress bar is not displayed

            out,err = proc.communicate()
        

        #If errors, raise an exception
        if proc.returncode!=0:
            err = str(err,'utf-8')
            raise RuntimeError(err)

        out = str(out,'utf-8').rstrip('\n')
        return out

    current_time=time.strftime('%d-%h--%H:%M',time.localtime(time.time()))
    tau_crit=[]
    tau_crit_err=[]
    nabla=[]
    nabla_err=[]
    import glob
    dat_files=sorted(glob.glob("data_*_*.dat"))
    boundaries=[]

    for file in dat_files:
        boundary=file.split("_")[1]
        boundaries.append(int(boundary))
    boundaries=sorted(boundaries)
    repeats=None 
    threads=None
    print(dat_files)
    #for file in dat_files:
    for boundary in boundaries: 
        file=list(glob.glob(f"data_{boundary}_*.dat"))[0]
        grep=system_call(f"grep -Po repeat=\\d+ {file}")
        if repeats and repeats!=int(grep.split("=")[1]):
            print("Warning: repeats differ between files, check the files. If this is intentional this can be ignored\n\
           This was noticed for {file}")

        repeats=int(grep.split("=")[1])

        grep=system_call(f"grep -Po threads=\\d+ {file}")
        if threads and threads!=int(grep.split("=")[1]):
            print(f"Warning: threads differ between files, check the files. If this is intentional this can be ignored\n\
            This was noticed for {file}")

        threads=int(grep.split("=")[1])
        


        boundary=file.split("_")[1]
        #boundaries.append(int(boundary)
        data=np.loadtxt(file)
        tau_crit.append(data[0][0]) #0 
        tau_crit_err.append(data[0][2]) #2
        
        nabla.append(data[0][1]) #0 
        nabla_err.append(data[0][3]) #2

    print(boundaries)
    def do_fit(boundaries,tau_crit,tau_crit_err):
        # print("tau error bars",tau_crit_err)
        try:
            par,cov = curve_fit(finitescaling,boundaries,tau_crit,sigma=tau_crit_err,absolute_sigma=True,maxfev=10000)
        except:
            return None
        b=np.linspace(boundaries[0],boundaries[-1],100)
        print("Val",boundaries[0],par[0],par[1],np.sqrt(np.diag(cov))[0],np.sqrt(np.diag(cov))[:])
        return par,cov
    def test(L,nu,a):
        return L**(-1/nu)*(a)
    def do_fit_nu(boundaries,nabla,nabla_err):
        # print("tau error bars",tau_crit_err)
        try:
            par,cov = curve_fit(test,boundaries,nabla,sigma=nabla_err,absolute_sigma=True,maxfev=10000)
        except:
            return None
        b=np.linspace(boundaries[0],boundaries[-1],100)
        print("Val_nu",boundaries[0],par[0],np.sqrt(np.diag(cov))[0],par,cov)#,np.sqrt(np.diag(cov))[:])
        return par,cov

    pars=[]
    bounds=[]
    for i in range(len(boundaries)):
        val1=do_fit_nu(boundaries[i:],nabla[i:],nabla_err[i:])
        if val1 is None:
            continue
        par,cov=val1

        nu=par[0]
        nu_err=np.sqrt(np.diag(cov))[0]
        #tau_crit_err_test=np.sqrt(np.array(tau_crit_err)**2)#+(boundaries**(-1/nu_err-1))**2)
        def test_fix_nu(L,phi_c,a):
            return L**(-1/nu)*(a)+phi_c
        def do_fit2(boundaries,tau_crit,tau_crit_err):
            from odrpack import odr_fit
            x=np.array(boundaries)**(-1/nu)
            x_err=1/nu**2*np.array(boundaries)**(-1/nu)*np.log(np.array(boundaries))*nu_err
            initial_guess=0.2895
            beta0=[10,0.2895]
            bounds=([0,0.2],[10,0.4])
            def linfit(x,beta):
                b1,b2=beta
                return b1*x+b2 
            fit=odr_fit(linfit,x,tau_crit,beta0,bounds=bounds,weight_x=1/x_err**2,weight_y=1/np.array(tau_crit_err)**2)
            print(fit.stopreason,fit.beta,fit.sd_beta, fit.beta[1]-fit.sd_beta[1],fit.beta[1]+fit.sd_beta[1]) 
            # print("tau error bars",tau_crit_err)
            # try:
            #     par,cov = curve_fit(test_fix_nu,boundaries,tau_crit,sigma=tau_crit_err,absolute_sigma=True,maxfev=10000)
            # except Exception as e:
            #     print(e)
            #     return None
            # b=np.linspace(boundaries[0],boundaries[-1],100)
            # print("Val",boundaries[0],par[0],par[1],np.sqrt(np.diag(cov))[0],np.sqrt(np.diag(cov))[:])
            

            return fit.beta,fit.sd_beta
        val=do_fit2(boundaries[i:],tau_crit[i:],tau_crit_err[i:])
        if val is not None:
            bounds.append(boundaries[i])
            #pars.append([par[0],np.sqrt(np.diag(cov))[0]])
            par,cov=val
            pars.append([par[1],cov[1]])

            #b=np.linspace(boundaries[0],boundaries[-1],100)
            #plt.scatter(boundaries[i:],tau_crit[i:],marker="o",s=5)
            #plt.errorbar(boundaries[i:],tau_crit[i:],yerr=tau_crit_err[i:],linestyle='')
            #plt.plot(b,test_fix_nu(b,*par),label=f"Fit,err={par[0],np.sqrt(np.diag(cov))[0]}")
            #ax=plt.gca()
            #ax.set_xlabel("L")
            #ax.set_ylabel(r"$\phi_c(L)$")
            #plt.legend()
            #plt.savefig(fname=f"fig_2_{current_time}_{ttype}.png")
            #exit()

        print("-------")

    print(pars)
    padrs=np.array(pars)
    plt.errorbar(bounds,padrs[:,0],padrs[:,1])
    plt.plot(np.linspace(bounds[0],bounds[-1],10),[0.2895]*10,linestyle="--",label="phi=0.2895")
    ax=plt.gca()
    plt.legend()
    plt.title("phi_c free nu")
    ax.set_ylim(0.288,0.293)
    ax.set_xlabel(r"L_min")
    ax.set_ylabel(r"phi_c")
    plt.savefig(fname="phi_c_driftin.png")
    plt.clf()
    # exit()
    comment = (f"Additional information:\n\
    total_repetitions={threads*repeats}\n\
    Data is stored in the following manner\n\
    critical value,error, boundaries\n")

if do_plot:
    b=np.linspace(boundaries[0],boundaries[-1],100)
    val1=do_fit_nu(boundaries[:],nabla[:],nabla_err[:])
    if val1 is None:
        raise ValueError("oof")
    parnu,covnu=val1

    nu=parnu[0]
    nu_err=np.sqrt(np.diag(covnu))[0]
    #tau_crit_err_test=np.sqrt(np.array(tau_crit_err)**2)#+(boundaries**(-1/nu_err-1))**2)
    def test_fix_nu(L,phi_c,a):
        return L**(-1/nu)*(a)+phi_c

    par,cov=do_fit2(boundaries,tau_crit,tau_crit_err)
    print(par,cov,"uhh")
    #plt.errorbar(boundaries,tau_crit,yerr=tau_crit_err,linestyle=None)
    # plt.plot(b,finitescaling(b,*par),label=f"Fit,err={par[0],np.sqrt(np.diag(cov))[0]}")
    b=b**(-1/nu)

    plt.errorbar(boundaries**(-1/nu),tau_crit,yerr=tau_crit_err,linestyle=None)
    plt.plot(b,par[0]*b+par[1],label=f"Fit,err={par[1],cov[1]}")
# np.savetxt(f"test.dat",np.array([[par[0],np.sqrt(np.diag(cov))[0],*boundaries]]),header=f"{comment}")
    ax=plt.gca()
    ax.set_xlabel("L")
    ax.set_ylabel(r"$\phi_c(L)$")

    plt.legend()

    plt.savefig(fname=f"fig_2_{current_time}_{ttype}.png")

if not read_data:
    os.system("rm -r ./tmp/")
