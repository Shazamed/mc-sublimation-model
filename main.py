import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import scipy

L = 8 # lattice size
kT = 1

lattice_array = np.ones([L,L,L])
surface_array = np.zeros([L,L,L])
surface_array[0,:,:] = surface_array[-1,:,:] = surface_array[:,0,:] = surface_array[:,-1,:] = surface_array[:,:,0] = surface_array[:,:,-1] = 1
lattice_list = []
plane_dict = {0:'[100]', 1:'[010]', 2:'[001]', 
              3:'[111]', 4:'[1-11]', 5:'[11-1]', 
              6:'[1-1-1]', 7:'[110]', 8:'[1-10]',
              9:'[101]', 10: '[10-1]', 11:'[011]',
              12:'[01-1]'} # 0 corresponds to x, 1 to y, 2 to z
# plane_lengths = {'[100]':[], '[010]':[], '[001]':[], '[111]':[], '[1-11]':[], '[11-1]':[], '[1-1-1]':[]}
plane_lengths = {}
for key in plane_dict:
    plane_lengths[plane_dict[key]] = []

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
    surface_array[particle_coord] = 0

def length_calc():
    non_zero_coords = np.where(lattice_array>0)
    for axis_idx in range(3):
        try:
            length = np.max(non_zero_coords[axis_idx])-np.min(non_zero_coords[axis_idx])+1
        except:
            length = 0
        plane_lengths[plane_dict[axis_idx]].append(length)
    
    # obtain long diagonals
    diagonal_basis = np.array([[1,1,1],[1,-1,1],[1,1,-1]])/(3**0.5)
    new_coords= np.matmul(diagonal_basis,non_zero_coords)
    for axis_idx in range(3,6):
        try:
            length = np.max(new_coords[axis_idx-3])-np.min(new_coords[axis_idx-3])
        except:
            length = 0
        plane_lengths[plane_dict[axis_idx]].append(length)
    
    axis_idx = 6
    diagonal_basis = np.array([[1,-1,-1],[1,-1,1],[1,1,-1]])/(3**0.5)
    new_coords= np.matmul(diagonal_basis,non_zero_coords)
    try:
        length = np.max(new_coords[0])-np.min(new_coords[0])
    except:
        length = 0
    plane_lengths[plane_dict[axis_idx]].append(length)

    # face diagonals
    axis_idx = 7
    face_bases_list = [
        np.array([[1,1,0],[1,-1,0],[0,0,1]])/(2**0.5),
        np.array([[1,0,1],[1,0,-1],[0,1,0]])/(2**0.5),
        np.array([[0,1,1],[0,1,-1],[1,0,0]])/(2**0.5)]
    for diagonal_basis in face_bases_list:
        new_coords= np.matmul(diagonal_basis,non_zero_coords)
        for i in range(2):
            try:
                length = np.max(new_coords[i])-np.min(new_coords[i])
            except:
                length = 0
            plane_lengths[plane_dict[axis_idx]].append(length)
            axis_idx += 1
        




i = 0
while True:
    i += 1
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
        lattice_list.append({"coords":np.copy(lattice_array),"iteration":i})
    length_calc()


# animation and plotting
for plane in plane_lengths.keys():
    plt.plot(plane_lengths[plane], label=plane)
# plt.plot(x_length, label="x")
# plt.plot(y_length, label="y")
# plt.plot(z_length, label="z")
plt.legend()
plt.title("Changes in length against number of iterations")
plt.show()


fig = plt.figure()

def updatefig(i):
    fig.clear()
    ax = fig.add_subplot(projection='3d')
    # ax.title("Monte Carlo Simulation for the Sublimation of a Crystal")
    x,y,z = lattice_list[i]["coords"].nonzero()
    ax.scatter(x,y,z,s=150)
    ax.set_xlim(0,L)
    ax.set_ylim(0,L)
    ax.set_zlim(0,L)
    plt.gcf().text(0.02,0.5,lattice_list[i]["iteration"],fontsize=14)
    plt.draw()
# quit()

print("writing animation")
writergif = animation.PillowWriter(fps=30)
anim = animation.FuncAnimation(fig, updatefig, frames=len(lattice_list),interval=200)
anim.save('simulation.gif', writer=writergif)
print("simulation.gif saved")
