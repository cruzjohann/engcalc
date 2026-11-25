"""Check against Chaudhry (2014) Appendix D: frequency response of a series piping system. Johann Cruz."""
import numpy as np
g = 9.81                                   # the FORTRAN uses 7.7047 = g*pi/4 with g = 9.81
tau0, H0, Q0, k, Tth = 1.0, 30.48, 0.0089, 0.2, 3.0
pipes = [(609.5, 0.61, 1219.0), (228.6, 0.30, 914.4)]     # L, D, a
wth = 2*np.pi/Tth
def response(wr):
    w = wr*wth; U = np.eye(2, dtype=complex)
    for L, D, a in pipes:                  # state {q, h}
        C = a/(g*np.pi*D**2/4); b = L/a
        F = np.array([[np.cos(b*w), -1j*np.sin(b*w)/C], [-1j*C*np.sin(b*w), np.cos(b*w)]])
        U = F @ U
    q1 = -(2*H0*k/tau0)/(U[1,0] - 2*H0/Q0*U[0,0])
    h, q = U[1,0]*q1, U[0,0]*q1
    return abs(h)/H0, abs(q)/Q0, np.degrees(np.angle(h)), np.degrees(np.angle(q))
book = [(.5,.037,.199,-95.257,-5.257),(1,.123,.190,-107.887,-17.887),(1.5,.076,.196,100.897,10.897),
        (2,.046,.199,-96.552,-6.552),(2.5,.149,.186,-111.940,-21.940),(3,.400,.000,179.998,89.998),
        (3.5,.149,.186,111.940,21.940),(4,.046,.199,96.551,6.551),(4.5,.076,.196,-100.898,-10.898),
        (5,.123,.190,107.886,17.886),(5.5,.037,.199,95.257,5.257),(6,.000,.200,-90.000,.000),
        (6.5,.037,.199,-95.258,-5.258),(7,.123,.190,-107.888,-17.888),(7.5,.076,.196,100.896,10.896),
        (8,.046,.199,-96.552,-6.552),(8.5,.149,.186,-111.941,-21.941),(9,.400,.000,179.995,89.995),
        (9.5,.149,.186,111.939,21.939),(10,.046,.199,96.551,6.551)]
print(f"{'wr':>5} {'H/H0 book':>9} {'mine':>7} {'Q/Q0 book':>9} {'mine':>7} {'phH book':>9} {'mine':>9}")
worst = 0
for wr, hb, qb, phb, pqb in book:
    h, q, ph, pq = response(wr)
    dph = (ph - phb + 180) % 360 - 180
    if hb > 0.001: worst = max(worst, abs(dph))
    print(f"{wr:5.1f} {hb:9.3f} {h:7.3f} {qb:9.3f} {q:7.3f} {phb:9.3f} {ph:9.3f}")
print('largest phase difference where |h| > 0:', round(worst, 3), 'deg')
