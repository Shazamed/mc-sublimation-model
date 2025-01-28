import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

L = 8 # lattice size
kT = 1

lattice_array = np.ones([L,L,L])
surface_array = np.zeros([L,L,L])
surface_array[0,:,:] = surface_array[-1,:,:] = surface_array[:,0,:] = surface_array[:,-1,:] = surface_array[:,:,0] = surface_array[:,:,-1] = 1
lattice_list = []

def calc_site_e(particle_coord):
    e_site = 0
    for axis in range(3):
        if particle_coord[axis] < L-1:
            neighbour_particle = list(particle_coord)
            neighbour_particle[axis] += 1
            if lattice_array[tuple(neighbour_particle)] == 1:
                e_site -= 1
        if particle_coord[axis] > 0:
            neighbour_particle = list(particle_coord)
            neighbour_particle[axis] -= 1
            if lattice_array[tuple(neighbour_particle)] == 1:
                e_site -= 1
    return e_site

def new_surface(particle_coord):
    for axis in range(3):
        if particle_coord[axis] < L-1:
            neighbour_particle = list(particle_coord)
            neighbour_particle[axis] += 1
            if surface_array[tuple(neighbour_particle)] != 1 and lattice_array[tuple(neighbour_particle)] != 0:
                surface_array[tuple(neighbour_particle)] = 1
        if particle_coord[axis] > 0:
            neighbour_particle = list(particle_coord)
            neighbour_particle[axis] -= 1
            if surface_array[tuple(neighbour_particle)] != 1 and lattice_array[tuple(neighbour_particle)] != 0:
                surface_array[tuple(neighbour_particle)] = 1
    surface_array[particle_coord] = 0 # value of 2 for removed particle


while True:
    e_site = 0
    try:
        random_particle = np.unravel_index(np.random.choice(np.flatnonzero(surface_array)), surface_array.shape)
    except:
        print("iterations complete")
        break
    e_site = calc_site_e(random_particle)
    probability = np.exp(e_site/kT)
    random_roll = np.random.random()
    if probability > random_roll:
        lattice_array[random_particle] = 0
        new_surface(random_particle)
        lattice_list.append(np.copy(lattice_array))

# animation and plotting
fig = plt.figure()

def updatefig(i):
    fig.clear()
    ax = fig.add_subplot(projection='3d')
    z,x,y = lattice_list[i].nonzero()
    ax.scatter(x,y,z,s=150)
    ax.set_xlim(0,L)
    ax.set_ylim(0,L)
    ax.set_zlim(0,L)
    plt.draw()


writergif = animation.PillowWriter(fps=30)
anim = animation.FuncAnimation(fig, updatefig, frames=len(lattice_list),interval=200)
anim.save('simulation.gif', writer=writergif)
