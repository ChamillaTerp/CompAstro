import numpy as np

# Solve for the HLL flux give to HD state vectors (left and right of the interface)
# stored in face[vars, 2, n,n,n], where 2=left, right
# face[:,] = [Density, P, vO, v1, v2], where vO is normal to the interface between
# left and right state, and v1, v2 are parallel 
def HLL(face,u):
    # shorthand values
    Dl = face[0,0]; Pl = face[1,0]; vl = face[2:5,0]
    Dr = face[0,1]; Pr = face[1,1]; vr = face[2:5,1]

    # Compute conserved variables [density, total energy, momentum]
    Etotl = Pl / (u.gamma - 1) + 0.5*Dl*np.sum(vl**2,axis=0)
    Cl = np.array([Dl,Etotl,Dl*vl[0],Dl*vl[1],Dl*vl[2]])

    Etotr = Pr / (u.gamma - 1) + 0.5*Dr*np.sum(vr**2,axis=0)
    Cr = np.array([Dr,Etotr,Dr*vr[0],Dr*vr[1],Dr*vr[2]])
    
    # Sound speed for each side of interface (l==left, r==right)
    c2_left  = u.gamma*Pl/Dl
    c2_right = u.gamma*Pr/Dr
    c_max = np.sqrt(np.maximum(c2_left,c2_right))

    # maximum wave speeds to the left and right (guaranteed to have right sign)
    SL = np.minimum(np.minimum(vl[0],vr[0])-c_max,0) # <= 0.
    SR = np.maximum(np.maximum(vl[0],vr[0])+c_max,0) # >= 0.

    # 1D HD fluxes for density, total energy, and momentum with advection velocity vl[0] / vr[0]
    #              Density    Etot             Momentum
    Fl = np.array([Dl*vl[0], (Etotl+Pl)*vl[0], Cl[2]*vl[0]+Pl, Cl[3]*vl[0], Cl[4]*vl[0]])
    Fr = np.array([Dr*vr[0], (Etotr+Pr)*vr[0], Cr[2]*vr[0]+Pr, Cr[3]*vr[0], Cr[4]*vr[0]])

    # HLL flux based on wavespeeds. If SL < 0 and SR > 0 (sub-sonic) then mix appropriately. 
    # Because conserved variables and fluxes are arrays we can simply loop
    Flux=np.empty_like(Fl) # shape (nv,n[0],n[1],n[2])
    for iv in range(u.nv):
        Flux[iv] = (SR * Fl[iv] - SL*Fr[iv] + SL*SR*(Cr[iv] - Cl[iv])) / (SR - SL)

    return Flux