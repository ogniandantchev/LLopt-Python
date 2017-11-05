import numpy as np
import matplotlib.pyplot as plt
import scipy.special as special
from scipy.integrate import quad
#from scipy import interpolate

V = 10.
n = 3.25 # s^-1
z = np.float64(2.)
D = 7.
#T = ..
#read*, KT
#KT= 0.33
KT1= 0.
rho= 1.2
#ni= 15e-6 # kinematic viscosity
R = D/2.
r1 = 1000000.
rns = 500000.
e1= 0.0001
J= V/(n*D)
P=13500 # 13500 Watts is close to the Betz limit for D = 7m, V = 10 m/s
cp = 8.*P / (np.pi * rho * D**2 * V**3)
# see https://en.wikipedia.org/wiki/Wind-turbine_aerodynamics or [4]
# for the definitions of the coefficients
tmp = np.sqrt(4.*cp/3.)
ct = -2.* tmp * np.sinh(np.arcsinh((-2.*cp**2)/tmp**3)/3.)  # check Ct sign!

#ct = 0.8888

KT = np.pi*ct*J**2/8.

print(f"Cp = {cp}, Ct = {ct}, KT = {KT}.")
print("Please note: Cp max is 16/27, or around 59%")

drh = 0.2 # dimensionless hub radius
rh= drh*R # hub radius, m
dh= rh * 2.  # hub diameter, m

omega = 2. * np.pi * n

# Chord distribution along the blade in 4 steps (DNV GL BLADED has 10)
chh= 0.08 #!14
ch1= 0.09 #!16
ch2= 0.07 #	!16
cht= 0.057 #	!14
chr1= 0.3
chr2= 0.7

deltah= 0.13
deltat= 0.04

B1=0.0733
B2=0.0174
# These are specific to the NACA 66 series mod a=0.8

xs= 0.
ts= 0.

# case where Circulation distribution gradient (dG/dx) is zero both at hub and tip of the blade
alfa= 0.5
beta= 0.5
# case where Circulation distribution gradient (dG/dx) is zero both at hub and tip of the blade
#alfa= -0.5
#beta= 0.5

lt = J/np.pi
#lt= V/(omega*R)
li= 0.944*lt-0.065 #!???
#li = .26*lt**2 + .812*lt + .35*KT + .03

xnaca = np.array([0., .0125, .025, .05, .075, .1, .2, .3, .4, .45, .5, .6, .7, .8, .9, .95, 1.])
fcnacaa08b005 = np.array([0., .0686, .142, .282, .389, .475, .725, .881, .97, .992, 1., .971, .877, .69, .352, .168, 0.])
ftnaca66 = np.array([0., .231, .306, .419, .508, .584, .8, .927, .99, 1., .992, .931, .807, .622, .375, .229, .006])
#integer(4), parameter:: maxft= 10

#np= 40
M= 21
No= M+1

#real(8) function chord(r)
def chord(r):
    #real(8) chh, ch1, ch2, cht, drh, chr1, chr2
    #common /cho/ chh, ch1, ch2, cht, drh, chr1, chr2
    if (r > chr2):
        return (r-chr2)*((cht-ch2)/(1.-chr2)) + ch2
    elif (r > chr1):
        return (r-chr1)*((ch2-ch1)/(chr2-chr1)) + ch1
    else:
        return (r-drh)*((ch1-chh)/(chr1-drh)) + chh
#end function

# x distance to dr
#real(8) function x2dr(x, drh)
#real(8) x, drh
def x2dr(x, drh):
    return (x*(1.-drh)+1.+drh)/2.
#end function

print("Omega= ", omega)

xkk = np.zeros((3, No)) #, order='F'
xii = np.zeros((3, No))

delta = np.zeros(No) #Chord distribution along the blade
c = np.zeros(No)

for k in range(0, No):
    xkk[0, k]= np.cos(np.pi*(2.*(k+1)-1.)/(2.*No+1.))
    xkk[1, k]= x2dr(xkk[0, k], drh)
    xkk[2, k]= R*xkk[1, k]
    #print(xkk[0,k], xkk[1,k], xkk[2,k])

print(xkk, "\n")

for i in range(0, No):
    xii[0, i]= np.cos(2*np.pi*(i+1)/(2*No+1.))
    xii[1, i]= x2dr(xii[0, i], drh)
    xii[2, i]= R*xii[1, i]
    #print(xii[0,k], xii[1,k], xii[2,k])
    
print(xii)



for i in range(0, No):
    delta[i]= (1.-np.log10(9.*(xii[1, i]-drh)/(1.-drh)+1.))* (deltah-deltat)+deltat
    c[i]= chord(xii[ 1, i])

print(delta, "\n", c)

plt.plot(xii[1,:], c)
plt.plot(xii[1,:], delta)
plt.show()

def itp(tan_beta_prime, r, rp):
    # see [4] and [1]
    s = np.sign(r - rp)
    if (s == 0):
        return np.sin(np.arctan(tan_beta_prime))
    else:
        rrr= rp*tan_beta_prime
        mu = r / rrr
        mup = rp / rrr
        #muh = rh / rrr
        u = np.sqrt(1. + mu**2)
        up = np.sqrt(1. + mup**2)
        #uh = np.sqrt(1. + muh**2)
     
        rer = 1. - rp / r
        ni = ((((u - 1.) / mu) * (mup / (up - 1.)))**(z * s)) * np.exp(z * abs(u - up))
        f1 = (s / (2. * z * mup)) * np.sqrt(up / u)
        f2 = 1. / (ni - 1.) - (s / (24. * z)) * ((9. * mup**2 + 2.) / up**3 + (3. * mu**2 - 2.) / u**3) * np.log(1. + (1. / (ni - 1.)))

        if (r > rp):
            return z * rer * (1. + 2 * mup * z * f1 * f2) 
        else:
            return 2 * (z**2) * rer * mup * f1 * f2 
    


def iap(tbp, r, rp):
    mu = r / (rp*tbp)
    return mu * itp(tbp, r, rp) + z * (rp / r - 1.)


#def yfun(u):
#    return u - .5 * np.log((u + 1.) / (u - 1.))


def Gfx(x):
    sum= 0.
    for i in range( 1, M+1):
        sum= sum + a[i] * special.eval_jacobi(i-1, alfa+1., beta+1., x) / np.float64(i)
    h= -a[0]*(.5*(beta-alfa+1.) * np.arccos(x)+(beta-alfa) * np.sqrt(1.-x**2.) + .5*(beta-alfa-1.) * x * np.sqrt(1-x**2.))
    # page 144 of [1]
    return h - .5*((1.-x)**(alfa+1.))*((1.+x)**(beta+1.))*sum


def Gdx(x):
    sum= 0.
    for i in range(0, M+1):
        sum= sum + a[i] * special.eval_jacobi(i, alfa, beta, x)

    return ((1.-x)**(alfa))*((1.+x)**(beta)) * sum


def ktpif(dr):
        #return interpolate.splev(dr, CSCOEF, der=0)
        
        be = np.arctan(lt/dr)
        bei= np.arctan(li/dr)
        u1= np.sin(be-bei) * np.sin(bei) / np.sin(be)
        #v1= 1 - np.tan(bei)*(np.pi*xii[1,i] / J + u1)
        #xdata[i]= xii[0,i]
        #tempg= Gfx(xdata[i])
        #tempvrm= np.pi*xii[1,i]/J + ((lt-li)/lt)*(xii[1,i]*li/(xii[1,i]**2+li**2))
        #tempvrm= np.pi*xii[1,i]/J + np.sin(be-bei) * np.sin(bei) / np.sin(be)
        cbi= np.cos(bei)
        
        vrstar = np.float64((np.pi*dr/J + u1) / cbi)
        
        return Gfx(dr) * vrstar * cbi * (1.-(Cdrag_fun(dr)/Clift_fun(dr)) * li / dr)

 
    
def dimdelta(x):
    return (1.-np.log10(9.*(x-drh)/(1.-drh)+1.))* (deltah-deltat)+deltat

def Clift_fun(dr):
    be = np.arctan(lt/dr)
    bei= np.arctan(li/dr)
    u1= np.sin(be-bei) * np.sin(bei) / np.sin(be)
    cbi= np.cos(bei)
    vrstar = np.float64((np.pi*dr/J + u1) / cbi)
    return 2.*np.pi*Gfx(dr)/(chord(dr)*vrstar)
        
    
def Cdrag_fun(dr):
    return .05808*(1. + 2.3*dimdelta(dr))/rns**(.1458)
    #rns= vstar * chord(dr*R) * R / ni # WTF защо Рейнолдс всеки път?
    #dimdelta = (1.-np.log10(9.*(dr-drh)/(1.-drh)+1.))* (deltah-deltat)+deltat
    #cd= .05808*(1. + 2.3*delta[i])/rns**(.1458)
    # for NACA 66 a=0.8 (mod b=0.05) airfoils only
    #cl[i]= 2.*np.pi*tempg/(c[i]*vri)

a = np.zeros( No )
s = np.zeros( (No, No) )
RHS = np.zeros(No)

phi = np.zeros(No)

#cl cd deltac
#allocate(s(1:M+1,0:M), a(0:M), RHS(1:M+1))
#allocate( xii(0:2,1:No), xkk(0:2,1:No))
#allocate( cl(1:No), cd(1:No), c(1:No))
#allocate( deltac(1:No), delta(1:No), phi(1:No))

#real(8)	xdata(8), fdata(8)
xdata = np.zeros(No)
fdata = np.zeros(No)
cl = np.zeros(No)
cd = np.zeros(No)
deltac = np.zeros(No)

CSCOEF = np.zeros(No)
#real(8)	break(8), cscoef(4,8)

# main loop for the Integro-differential equation
KT = np.pi*ct*J**2/8.
KT1 = 0.

while (abs(KT - KT1) > e1):
    for i in range(0, No):
        for nn in range( 0, No):
            sum= 0.
            for k in range(0, No):
                c2tk2= np.cos(np.pi*(2.*k-1.)/(2.*(2.*No+1.)))**2
                tif= itp(li/xkk[1,k], xii[1,i], xkk[1,k])
                sum= sum + c2tk2*special.eval_jacobi(nn, -.5, .5, xkk[0,k])*tif/(xkk[0,k]-xii[0,i])
                # from SciPy -- special.eval_jacobi(nn, alpha, beta, x)
                #end do
            s[i,nn]= (1/(1.-drh))*(4.*np.pi/(2.*No+1))*sum
            #end do
        s[i,0]= s[i,0] + np.pi*z/(2.*xii[1,i])
        
        RHS[i]= 2.*((lt-li)/lt)*(xii[1,i]*li/(xii[1,i]**2+li**2))
        #be = np.arctan(lt/xii[1,i])
        #bei= np.arctan(li/xii[1,i])
        #RHS[i]= np.sin(be-bei) * np.sin(bei) / np.sin(be)
        #end do

    #print(s, "\n",RHS)
    #call dlsarg(M+1, s, M+1, RHS, 1, a)
    #Solves a real general system of linear equations with iterative refinement (IMSL again).
    # s*a=RHS ... solve(s,RHS)
    a = np.linalg.solve(s, RHS)
    
    RES1 = quad(ktpif, drh, 1.)

#    for i in range( 0, No):
#        be = np.arctan(lt/xii[1,i])
#        bei= np.arctan(li/xii[1,i])
#        u1= np.sin(be-bei) * np.sin(bei) / np.sin(be)
#        v1= 1 - np.tan(bei)*(np.pi*xii[1,i] / J + u1)
#        xdata[i]= xii[0,i]
#        tempg= Gfx(xii[0,i])
#        #tempg= Gfx(xdata[i])
#        tempvrm= np.pi*xii[1,i]/J + ((lt-li)/lt)*(xii[1,i]*li/(xii[1,i]**2+li**2))
#        #tempvrm= np.pi*xii[1,i]/J + np.sin(be-bei) * np.sin(bei) / np.sin(be)
#        cbi= np.cos(bei)
#        
#        vri= np.sqrt((1-v1)**2 + (np.pi*xii[1,i]/J)**2)
#        rns= vri * c[i] * R / ni # WTF защо Рейнолдс всеки път?
#        #rns= 500000
#
#        cd[i]= .05808*(1. + 2.3*delta[i])/rns**(.1458)
#        # for NACA 66 a=0.8 (mod b=0.05) airfoils only
#
#        cl[i]= 2.*np.pi*tempg/(c[i]*vri)
#        
#        fdata[i]= tempg*vri*cbi*(1.-(cd[i]/cl[i]) * li / xii[1,i])
#        #fdata[i]= tempg*vri*cbi*(1.-(cd[i]/cl[i]) * li / xii[1,i])
#        #end do

    #CALL dCSINT (No, XDATA, FDATA, BREAK, CSCOEF)
    #This IMSL function computes the cubic spline interpolant with the ‘not-a-knot' condition.

    #CSCOEF = interpolate.splrep(xdata, fdata, s=0)
    
    #p.interpolate.LSQUnivariateSpline(x, y, t, w=None, bbox=[None, None], k=3, ext=0, check_finite=False)

    #CALL dQDAGS (ktpif, drh, 1., 0., 0.0001, RES1, ERR1)
    # from scipy.integrate import quad
    
    #if2 = interpolate.InterpolatedUnivariateSpline(xdata, fdata, w=None, bbox=[None, None], k=3, ext=0, check_finite=False)
    #RES1 = quad(if2, drh, 1.)
    
    #RES1 = quad(ktpif, drh, 1.)
    

    #print(RES1)
    #print("ERR= ", RES1[1])
    
    #RES1 = quad(ktpif, drh, 1.)

    #KT1 = -z * J**2. * np.pi * RES1[0] / 2.
    KT1 = np.pi * z * J**2 * RES1[0] / 2.

    print( KT1, li)
    print(f"KT = {KT}, KT1 ={KT1}, ERR ={KT-KT1}")

    li= li * KT / KT1
#end do

#CALL dWRRRN ('a', 1, M+1, a, 1, 0)
# from IMSL again - now it's simply print...
print(f"a = { a}, lambda_t ={lt}, lambda_i ={li}")

for i in range(0, No):
    deltac[i]= B1 * Clift_fun(xii[1,i])
    phi[i]= np.arctan(li / xii[1,i]) - B2 * Clift_fun(xii[1,i])
    # good candidate for future rework -- need to read real airfoilds from data files
#        dimdelta = (1.-np.log10(9.*(dr-drh)/(1.-drh)+1.))* (deltah-deltat)+deltat

print("Phi = ", phi)
plt.figure()
plt.plot(xii[1,:], phi)
plt.show()

print("DeltaC =", deltac)
plt.figure()
plt.plot(xii[1,:], deltac)
plt.show()


#plt.figure()
#plt.scatter(x, Clift_fun(x))
#plt.scatter(x, Cdrag_fun(x))
#plt.show()

#plt.figure()
#plt.plot(xdata, fdata)
#plt.show()

def two_scales(ax1, x, data1, data2, c1, c2):
    """
    Parameters
    ----------
    ax : axis - Axis to put two scales on
    x : array-like - x-axis values for both datasets
    data1: array-like - Data for left hand scale
    data2 : array-like - Data for right hand scale
    c1, c2 : color - Colors for line 1 and 2
    

    Returns
    -------
    ax : axis
        Original axis
    ax2 : axis
        New twin axis
    """
    ax2 = ax1.twinx()

    ax1.plot(x, data1, color=c1)
    ax1.set_xlabel('r/R')
    ax1.set_ylabel('G')

    ax2.plot(x, data2, color=c2)
    ax2.set_ylabel('dG/dr')
    return ax1, ax2

# Change color of each axis
def color_y_axis(ax, color):
    """Color for the axes."""
    for t in ax.get_yticklabels():
        t.set_color(color)
    return None


x = np.linspace(drh, 1, 70)
yg = Gfx(x)
ygp = Gdx(x)

# Create axes
fig, ax = plt.subplots()
ax1, ax2 = two_scales(ax, x, yg, ygp, 'r', 'b')
color_y_axis(ax1, 'r')
color_y_axis(ax2, 'b')
plt.show()


#x = np.zeros(No)
#x = xii[0,:]

