# ---------------------------------------------------------------------------- #
# Import relevant Python libraries                                             #
# ---------------------------------------------------------------------------- #

import numpy as np
import os
import time

# ---------------------------------------------------------------------------- #
# Basic definitions                                                            #
# ---------------------------------------------------------------------------- #

# Make a snapshot directory

os.mkdir('SNAPS') 

# Define the number of particles

n_part = 50 

# Define the box length and dimension

l_box = 10.0
d_sys = 2

# Define the particle mass

mass = 1.0

# Define the thermal energy

kBT = 1.0

# Define the WCA parameters

sig = 1.0
eps = 1.0
cut = (2.0**(1.0/6.0))*sig
preforce = 24.0*eps/sig

# Equipartition prefactor

preequip = np.sqrt(d_sys*kBT/mass)

# Define the time step

dt = 0.001

# Define the number of steps 

n_step = 10000
n_out = 100
n_print = 10

# ---------------------------------------------------------------------------- #
# Initialize the positions, velocities, and acceleration                       #
# ---------------------------------------------------------------------------- #

# Set the seed (read in the run's PID and time for more randomness)

pid = int(os.getpid())
epoch = int(time.time())
pseed = int(pid*epoch%(2**32))
np.random.seed(pseed)

# Particle positions need to be randomized

posit = np.zeros([n_part,d_sys])
rij = np.zeros(d_sys)

for i in range(0,n_part):

    coords = np.zeros(d_sys)

    dist2 = sig*sig

    # remove overlaps

    while( dist2 < 1.05*cut*cut):

        for j in range(0,d_sys):

            coords[j] = l_box*np.random.rand()

        min_dist = l_box*l_box

        # Check the minimum distance with previous particles

        for j in range(0,i):

            # Calculate the distance vector

            for k in range( 0, d_sys ):
                
                rij[k] = posit[j,k] - coords[k]
                
                # deal with the periodic condition 
                
                while ( rij[k] > 0.5*l_box ):
                    
                    rij[k] -= l_box

                while ( rij[k] <= -0.5*l_box ):
                    
                    rij[k] += l_box

            # Compute length of the vector
                    
            dist2 = np.dot(rij,rij)

            if (dist2 < min_dist):
                min_dist = dist2

        dist2 = min_dist
        
    posit[i] = coords

# Set the velocities according to equipartition

veloc = np.zeros([n_part,d_sys])

# Random positions

for i in range(0,n_part):

    vel = np.zeros(d_sys)

    for j in range(0,d_sys):

        vel[j] = (np.random.rand() - 0.5)

    veloc[i] = vel

# Remove center of mass velocity

v_cm = np.zeros(d_sys)

for i in range(0,n_part):

    v_cm += veloc[i]

v_cm /= n_part

for i in range(0,n_part):

    veloc[i] -= v_cm

# Set the initial kinetic energy to equipartition

e_kin = 0.0

for i in range(0,n_part):

    e_kin += np.dot(veloc[i],veloc[i])

e_kin *= 0.5*mass
e_equi = (d_sys/2.0)*kBT*n_part

for i in range(0,n_part):

    veloc[i] *= np.sqrt(e_equi/e_kin)
    
e_kin = e_equi

# Compute the forces acting on the particles ----------------------------- #

e_pot = 0.0
p_vir = 0.0
force = np.zeros([n_part,d_sys])

# Loop over all particle pairs

for i in range( 0, n_part ):
    
    for j in range( i+1, n_part ):

        # Calculate the distance vector

        for k in range( 0, d_sys ):
            
            rij[k] = posit[i,k] - posit[j,k]
            
            # deal with the periodic condition 
            
            while ( rij[k] > 0.5*l_box ):
                
                rij[k] -= l_box

            while ( rij[k] <= -0.5*l_box ):
                
                rij[k] += l_box

        # Compute length of the vector
                
        dist = np.sqrt(np.dot(rij,rij))
        
        if ( dist < cut ):
            
            # Calculate the Weeks-Chandler-Anderson force

            sir1 = sig/dist
            sir3 = sir1*sir1*sir1
            sir6 = sir3*sir3
            prefactor = preforce*sir1*sir6*( 2.0*sir6 - 1.0 )
            preenergy = 4.0*eps*( sir6*sir6 - sir6 + 0.25 )
                
        else:
            
            prefactor = 0.0
            preenergy = 0.0
            
        rhat = (rij/dist)
        force[i] += prefactor*rhat
        force[j] -= prefactor*rhat

        # Add to the potential energy

        e_pot += preenergy

        # Add to the pressure

        p_vir += np.dot(force[i],rij)

# ---------------------------------------------------------------------------- #
# Loop through the simulation                                                  #
# ---------------------------------------------------------------------------- #

t0 = time.time()

for i_time in range( 0, n_step ):

    # Output the system's properties ----------------------------------------- #
            
    e_tot = e_kin + e_pot
    temp = np.sqrt(2.0*e_kin/(d_sys*n_part))

    if ( i_time % n_out == 0 ):
        
        vol = l_box**d_sys
        pres = kBT*n_part/vol + p_vir/(vol*d_sys)  
        
        print('time: {} out of {}, ekin: {}, epot: {}, etot: {}, temp: {}, pres: {}'.format(i_time, n_step, e_kin, e_pot, e_kin + e_pot, temp, pres))

    # Output the particle coordinates

    if ( i_time % n_print == 0 ):
        
        np.savetxt('SNAPS/pos_{}.dat'.format(i_time/n_print), posit, delimiter=' ')

    # The first Verlet half-time step ---------------------------------------- #

    # Calculate the new particle positions

    for i in range( 0, n_part ):
        
        for j in range( 0, d_sys ):
            
            posit[i,j] += veloc[i,j]*dt + 0.5*(force[i,j]/mass)*dt*dt
            
            # deal with the periodic condition 
                
            while ( posit[i,j] > l_box ):
                
                posit[i,j] -= l_box

            while ( posit[i,j] < 0.0 ):
                
                posit[i,j] += l_box

    # Calculate the new particle velocities

    for i in range( 0, n_part ):
        
        for j in range( 0, d_sys ):
            
            veloc[i,j] += 0.5*dt*(force[i,j]/mass) 

    # Compute the forces acting on the particles ----------------------------- #

    force = np.zeros([n_part,d_sys])
    e_pot = 0.0
    p_vir = 0.0

    # Loop over all particle pairs

    for i in range( 0, n_part ):
        
        for j in range( i+1, n_part ):

            # Calculate the distance vector

            for k in range( 0, d_sys ):
                
                rij[k] = posit[i,k] - posit[j,k]
                
                # deal with the periodic condition 
                
                while ( rij[k] > 0.5*l_box ):
                    
                    rij[k] -= l_box

                while ( rij[k] <= -0.5*l_box ):
                    
                    rij[k] += l_box

            # Compute length of the vector
                    
            dist = np.sqrt(np.dot(rij,rij))
            
            if ( dist < cut ):
                
                # Calculate the Weeks-Chandler-Anderson force
    
                sir1 = sig/dist
                sir3 = sir1*sir1*sir1
                sir6 = sir3*sir3
                prefactor = preforce*sir1*sir6*( 2.0*sir6 - 1.0 )
                preenergy = 4.0*eps*( sir6*sir6 - sir6 + 0.25 )
                    
            else:
                
                prefactor = 0.0
                preenergy = 0.0
                
            rhat = (rij/dist)
            force[i] += prefactor*rhat
            force[j] -= prefactor*rhat

            # Add to the potential energy

            e_pot += preenergy

            # Add to the pressure

            p_vir += np.dot(force[i],rij)

    # Second half step of the Velocity-Verlet algorithm ---------------------- #

    # Calculate the new particle velocities

    e_kin = 0.0
    for i in range( 0, n_part ):
        
        for j in range( 0, d_sys ):
            
            veloc[i,j] += 0.5*dt*(force[i,j]/mass)

        e_kin += np.dot(veloc[i],veloc[i])

    e_kin *= 0.5*mass

t1 = time.time()
del_t = t1 - t0
print("elapsed time: %10.3e" % del_t )
