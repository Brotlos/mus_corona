import numpy as np
import matplotlib.pyplot as plt


# Calculate Changes in S,I,R,X,D groups
def sirxd_update_population(S, I, R, Xs, Xi, Dn, Di, β, γ, δ, κs, κi, κe, v=0.0, μ=0.0, dt=1.0):
    """Calculate Changes in S,I,R,X,D groups

    Args:
        S (float): Number of susceptible individuals
        I (float): Number of infectious individuals
        R (float): Number of recovered individuals

        Xs (float): Number of susceptible individuals in quarantine
        Xi (float): Number of infectious individuals in quarantine
        Dn (float): Number of naturally deceased individuals
        Di (float): Number of deceased infectious individuals

        β (float): Infection rate
        γ (float): Recover rate
        δ (float): Death rate

        κs (float): Quarantine rate for susceptible individuals #TODO ω nach Ebola-Modell
        κi (float): Quarantine rate for infectious individuals #TODO ω nach Ebola-Modell
        κe (float): Quarantine-End rate for susceptible individuals #TODO ω nach Ebola-Modell

        v (float): (natural) Birth rate
        μ (float): (natural) Death rate

        dt (float): Duration of time step

    Returns:
        (float, float, float, float, float, float, float): 
        Returns updated size of each S,I,R,X,D group
    """



    # Population Size
    N = S + I + R + Xs + Xi

    # Change in susceptible individuals
    dS = (v*N +κe*Xs -β*S*I/N -κs*S -μ*S) * dt

    # Change in infectious individuals
    dI = (β*S*I/N -γ*I -κs*I -κi*I -δ*I -μ*I) * dt

    # Change in recovered individuals
    dR = (γ*I +γ*Xi -μ*R) * dt

    # Change in susceptible individuals in quarantine
    dXs = (κs*S -κe*Xs -μ*Xs) * dt

    # Change in infectious individuals in quarantine
    dXi = (κs*I +κi*I -γ*Xi -δ*Xi -μ*Xi) * dt

    # Change in deceased individuals (cause: natural death)
    dDn = (μ*S +μ*Xs +μ*I +μ*Xi +μ*R) * dt

    # Change in deceased infectious individuals (cause: epidemic)
    dDi = (δ*I +δ*Xi) * dt

    # Update Population
    S  += dS
    I  += dI
    R  += dR
    Xs += dXs
    Xi += dXi
    Dn += dDn
    Di += dDi

    # Lower Bound for group size
    if S < 0.0:
        S = 0.0
    if I < 0.0:
        I = 0.0
    if R < 0.0:
        R = 0.0
    if Xs < 0.0:
        Xs = 0.0
    if Xi < 0.0:
        Xi = 0.0
    if Dn < 0.0:
        Dn = 0.0
    if Di < 0.0:
        Di = 0.0
    N = S+I+R+Xs+Xi

    # Retun new S,I,R,X,D values
    return S, I, R, Xs, Xi, Dn, Di, N

def sirxd_update_groups(groups, rates, dt=1.0):
    return sirxd_update_population(*groups, *rates, dt)

# Simulate epidemic using SIRXD model and constant rates
def sim_epidemic_sirxd(N, I, T, rates, time_step=1.0, adapt_birthrate=True):

    # Init
    cS, cI, cR, cXs, cXi, cDn, cDi, cN = N-I, I, 0, 0, 0, 0, 0, 0
    S, I, R, Xs, Xi, Dn, Di, N = [cS], [cI], [cR], [cXs], [cXi], [cDn], [cDi], [cS+cI+cR]
    β, γ, δ, κs, κi, κe, v, μ = rates
    dt = time_step

    # Simulate
    time = list(range(1, T))
    for _ in time:

        if adapt_birthrate is True:
            # adapt birth rate to match effective death rate
            cN = cS+cI+cR+cXs+cXi
            v = μ + (δ*(cI+cXi)*dt)/cN 
        
        # update groups
        cS, cI, cR, cXs, cXi, cDn, cDi, cN = sirxd_update_population(
            cS, cI, cR, cXs, cXi, cDn, cDi, 
            β, γ, δ, κs, κi, κe, v, μ, 
            dt
        )

        # Save current time step
        S.append(cS)
        I.append(cI)
        R.append(cR)
        Xs.append(cXs)
        Xi.append(cXi)
        Dn.append(cDn)
        Di.append(cDi)
        N.append(cN)

    # Plot data
    time = [0] + time
    sirxd_plot(time, S, I, R, Xs, Xi, Dn, Di, N)

def sirxd_plot(time, S, I, R, Xs, Xi, Dn, Di, N):
    """Plots the given SIRXD model given by it's groups/classes"""
    plt.plot(
        time, S,
        time, I,
        time, R,
        time, Xs,
        time, Xi,
        time, Dn,
        time, Di,
        time, N
    )
    plt.xlabel('Time')
    plt.ylabel('S,I,R,X,D')
    plt.xlim((0, time[-1]))# plt.xlim((0, T))
    plt.ylim((0, 1.5*N[0]))
    plt.legend((
        "Susceptible", 
        "Infectious", 
        "Recovered", 
        "Susceptible in quarantine", 
        "Infectious in quarantine", 
        "Naturally deceased", 
        "Deceased infectious", 
        "Total population"))
    plt.show()



if __name__ == "__main__":

    # set rates  β,      γ,      δ,      κs,     κi,    κe,    v,       μ
    rates =     (0.5,    0.03,   0.01,   0.01,   0.1,   0.1,   0.0,    0.005)

    # Start simulation
    sim_epidemic_sirxd(N=1000, I=3, rates=rates, T=300)
