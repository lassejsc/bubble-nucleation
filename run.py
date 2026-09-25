
import numpy as np
import matplotlib.pyplot as plt
import os
import subprocess
import time
import scipy.stats as st
from scipy.optimize import curve_fit
from scipy.special import erf
import argparse

#This is the original code used in the thesis, I might update this later to be more clear, but I recommend making your own to use with the program
#----------------------------------------------------------------
#   Input parameters
#
#
#----------------------------------------------------------------


parser= argparse.ArgumentParser(add_help=False)
parser.add_argument("-t","--transition",help="The type of transition to use, current options: fast,slow,equal_size",required=True)


args,_ = parser.parse_known_args()
ttype=args.transition.lower()
if ttype=="fast":
    parser.add_argument("-G_f",help="value of percolation rate Gamma(t_f)=Gamma_f",required=True)
    parser.add_argument("-tau_f",help="maximum value of tau parameter, tau=beta*(t-t_f)",required=True)
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

parser.add_argument(
"-h", "--help",
action="help",
default=argparse.SUPPRESS,
help="show this help message and exit"
)
args=parser.parse_args()


#Main repeats
repeats=args.n_repeat
threads=args.threads

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


#boundary=30
#boundaries = np.arange(20,61,10)
boundaries=args.boundaries
print(boundaries)



#----------------------------------------------------------------
#   File manipulation, the data gets saved to /tmp/data.dat momentarily
#   if this folder already exists throws exception
#----------------------------------------------------------------

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
    return a*L**(-1/nu)+phi_c

#----------------------------------------------------------------
#   Main loop
#
#
#----------------------------------------------------------------

#For initial value forking
repeat_step=100
repeat_fork=4



if ttype=="equal_size":
    tau_div=12
    tau_upper=0.5
    tau_lower=0.2
    adjust_size=0.01
else:
    tau_div=6
#Initial starting points to refine from
    tau_upper=0
    tau_lower=-2
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
tau=(tau_upper+tau_lower)/2
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

    print(f"{(time.time()-start)/60}s")

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
       #error_max.append(error)
        error_both.append(error)
        #error_max.append(max(error))
        if probability>=0.35 and tau_l<tau:
            tau=tau_l
    data=np.array([taus,probabilities])
    data=data[:,np.argsort(data[0])]
    mask=np.logical_and(data[1]>=0.35,data[1]<=0.65) #Mask out data not in correct range
    taus=data[0][mask]
    probabilities=data[1][mask]
    error_max=np.array(error_max)[mask]
    error_both=np.array(error_both)[mask]
    print(probabilities)

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
plt.savefig(fname=f"fig_1_{current_time}_{simulation_name}.png")

plt.figure(2)

plt.errorbar(boundaries,tau_crit,yerr=tau_crit_err)
print(tau_crit_err)
par,cov = curve_fit(finitescaling,boundaries,tau_crit,sigma=tau_crit_err,absolute_sigma=True)
b=np.linspace(20,100,100)
plt.plot(b,finitescaling(b,*par),label=f"Fit,err={par[0],np.sqrt(np.diag(cov))[0]}")
ax=plt.gca()
ax.set_xlabel("L")
ax.set_ylabel(r"$\phi_c(L)$")

plt.legend()

plt.savefig(fname=f"fig_2_{current_time}_{simulation_name}.png")

os.system("rm -r ./tmp/")
