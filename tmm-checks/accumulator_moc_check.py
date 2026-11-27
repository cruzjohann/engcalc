"""
Time-domain check of accumulator_example_8_2.py by the method of characteristics.

Shares no code with the transfer matrix scripts. Example 8-2 (Chaudhry, Applied
Hydraulic Transients, 3rd ed.), frictionless, with:
  - the full valve equation Q = Cv*tau*sqrt(H), tau = 1 + k*sin(w t), discharging to atmosphere
  - a gas-charged accumulator at the tee or at the valve, gas law (H + Hatm)*V^n = const
The head amplitude at the forcing frequency is taken from the last part of the run.
Johann Cruz
"""

import numpy as np, sys
from scipy.optimize import brentq
g=9.80665; rho=999.0; patm=101325.0; Hatm=patm/(rho*g)
Ho,Qo,tau0=100.0,0.314,1.0
P={'1':(500.0,1000.0,1.81),'2':(250.0,1000.0,1.05),'b':(200.0,800.0,1.36)}
dt=0.0125
def mk(key):
    L,a,D=P[key]; A=np.pi*D*D/4; N=int(round(L/a/dt)); return dict(N=N,B=a/(g*A))
def run(f, k=0.2, V=0.0, where=None, n=1.0, T=420.0, Tavg=120.0):
    p1,p2,pb=mk('1'),mk('2'),mk('b')
    H1=np.full(p1['N']+1,Ho); Q1=np.full(p1['N']+1,Qo)
    H2=np.full(p2['N']+1,Ho); Q2=np.full(p2['N']+1,Qo)
    Hb=np.full(pb['N']+1,Ho); Qb=np.zeros(pb['N']+1)
    Cv=Qo/(tau0*np.sqrt(Ho)); w=2*np.pi*f
    Kgas=(Ho+Hatm)*V**n if V>0 else 0; Vg=V; Qa=0.0
    nsteps=int(T/dt); rec=[]; ts=[]
    for it in range(1,nsteps+1):
        t=it*dt
        tau=tau0+k*np.sin(w*t)
        # interior
        def interior(H,Q,B):
            Cp=H[:-2]+B*Q[:-2]; Cm=H[2:]-B*Q[2:]
            Hn=H.copy(); Qn=Q.copy()
            Hn[1:-1]=0.5*(Cp+Cm); Qn[1:-1]=(Cp-Cm)/(2*B); return Hn,Qn
        H1n,Q1n=interior(H1,Q1,p1['B']); H2n,Q2n=interior(H2,Q2,p2['B']); Hbn,Qbn=interior(Hb,Qb,pb['B'])
        # reservoir
        H1n[0]=Ho; Q1n[0]=(Ho-(H1[1]-p1['B']*Q1[1]))/p1['B']
        # dead end
        Qbn[-1]=0.0; Hbn[-1]=Hb[-2]+pb['B']*Qb[-2]
        # tee
        Cp1=H1[-2]+p1['B']*Q1[-2]; Cm2=H2[1]-p2['B']*Q2[1]; Cmb=Hb[1]-pb['B']*Qb[1]
        S=1/p1['B']+1/p2['B']+1/pb['B']; R=Cp1/p1['B']+Cm2/p2['B']+Cmb/pb['B']
        if where=='tee':
            def F(Hp):
                Qan=R-Hp*S
                Vn=Vg-0.5*dt*(Qa+Qan)
                return (Hp+Hatm)*Vn**n-Kgas
            Hp=brentq(F,H1[-1]-60,H1[-1]+60,xtol=1e-10); Qan=R-Hp*S; Vg=Vg-0.5*dt*(Qa+Qan); Qa=Qan
        else:
            Hp=R/S
        H1n[-1]=Hp; Q1n[-1]=(Cp1-Hp)/p1['B']; H2n[0]=Hp; Q2n[0]=(Hp-Cm2)/p2['B']; Hbn[0]=Hp; Qbn[0]=(Hp-Cmb)/pb['B']
        # valve
        Cp2=H2[-2]+p2['B']*Q2[-2]; B2=p2['B']
        if where=='valve':
            def F(Hp):
                Qin=(Cp2-Hp)/B2; Qan=Qin-Cv*tau*np.sqrt(max(Hp,0))
                Vn=Vg-0.5*dt*(Qa+Qan)
                return (Hp+Hatm)*Vn**n-Kgas
            Hp=brentq(F,max(H2[-1]-60,1e-6),H2[-1]+60,xtol=1e-10)
            Qin=(Cp2-Hp)/B2; Qan=Qin-Cv*tau*np.sqrt(Hp); Vg=Vg-0.5*dt*(Qa+Qan); Qa=Qan
        else:
            c=Cv*tau*B2; x=(-c+np.sqrt(c*c+4*Cp2))/2; Hp=x*x; Qin=(Cp2-Hp)/B2
        H2n[-1]=Hp; Q2n[-1]=Qin
        H1,Q1,H2,Q2,Hb,Qb=H1n,Q1n,H2n,Q2n,Hbn,Qbn
        if t>T-Tavg: rec.append(Hp); ts.append(t)
    ts=np.array(ts); h=np.array(rec)-Ho
    nper=int(Tavg*f); m=ts>=ts[-1]-nper/f+1e-9
    A=2*np.mean(h[m]*np.exp(-1j*w*ts[m]))
    return 2*abs(A)/Ho
if __name__=="__main__":
    print("Head response hr at the valve")
    cases=[("no accumulator",0.667,0,None,0.1688),("0.5 m3 at tee",0.667,0.5,'tee',0.187),("2 m3 at tee",0.667,2.0,'tee',0.211),
           ("0.5 m3 at valve",0.667,0.5,'valve',0.106),("2 m3 at valve",0.667,2.0,'valve',0.018),("0.5 m3 at valve, 1.0 Hz",1.0,0.5,'valve',0.044)]
    for name,f,V,wh,tmm in cases:
        r=run(f if f!=0.667 else 4.189/(2*np.pi),V=V,where=wh)
        print(f"{name:26s} TMM {tmm:.3f}  MOC {r:.3f}  diff {100*(r/tmm-1):+.1f}%",flush=True)
    print("\nShifted resonances, k = 0.02 (linear range)")
    for name,V,wh,fs in [("0.5 m3 at valve (TMM 0.541 Hz)",0.5,'valve',np.arange(0.533,0.5501,0.002)),
                         ("0.5 m3 at tee (TMM 0.282 Hz)",0.5,'tee',np.arange(0.276,0.2881,0.0015))]:
        r=[run(f,k=0.02,V=V,where=wh,T=900,Tavg=200) for f in fs]
        i=int(np.argmax(r)); print(f"{name:32s} MOC peak at {fs[i]:.3f} Hz, hr = {r[i]:.3f}",flush=True)
