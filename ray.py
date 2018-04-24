import matplotlib as mpl
mpl.use('Agg')

import yt
from yt import YTArray
import trident
import numpy as np
import numpy.random as npr
import matplotlib.pyplot as plt
import time
from mpi4py import MPI
import os
import shutil

from yt.funcs import mylog
mylog.setLevel(40)  # only logs errors and critical problems

comm = MPI.COMM_WORLD

t_start = time.time()

# number of pixels per dimension (for Aitoff projection)
pixels_per_dim = 512

# set the number of kpc you want to remove from the start of the ray
# trying to remove the data from the disk
remove_first_N_kpc = 0.0

# Load dataset
fn = "/mnt/scratch/dsilvia/simulations/reu_sims/MW_1638kpcBox_800pcCGM_200pcDisk_lowres/DD1500/DD1500"
ds = yt.load(fn)

# Add H I & O VI ion fields using Trident
trident.add_ion_fields(ds, ions=['O VI', 'H I'], ftype="gas")

# Specify position where the ray starts (position in the volume)
c = ds.arr([0.5, 0.5, 0.5], 'unitary')
c = c.in_units('kpc')

# location in the disk where we're setting our observations
# move start of ray to 8kpc from the center of the disk (about where the Sun is in the Milky Way)
X = YTArray([8., 0., 0.], 'kpc')
ray_start = c - X

# Length of Ray
R = 18.0

# do you want projections of the spheres?  True/False
MakeProjections = False

# do you want debug information while the calculation goes on?  True/False
Debug = True

# do you want to write out ray files?  True/False
# must be True if you don't already have ray files saved
WriteRays = True


# fraction of rays that make plots?  (determined randomly)
# plots of H I density, O VI density, density, temperature, metallicity, line of sight velocity
fraction_ray_plots = 0.0

# Create dictionaries for making projections
dens_dict = {'field':'density',  'file_name':'Density',  'title':'Density', 'colormap':'kelp'}
temp_dict = {'field':'temperature',  'file_name':'Temperature',  'title':'Temperature', 'colormap':'RdBu_r'}
metal_dict = HI_dict = {'field':'metallicity',  'file_name':'Metallicity',  'title':'Metallicity', 'colormap':'arbre'}
HI_dict = {'field':'H_p0_number_density',  'file_name':'HI_density',  'title':'HI Density', 'colormap':'plasma'}
OVI_dict = {'field':'O_p5_number_density',  'file_name':'OVI_density',  'title':'OVI Density', 'colormap':'viridis'}

field_dicts = [dens_dict, temp_dict, metal_dict, HI_dict, OVI_dict]

# sets the width of the projections being made
proj_scales = [40., 400.]

proj_dims = ['x','y','z']

# this creates a folder called ray_files in which the .h5 ray files will be saved, if you don't already have one
if comm.rank == 0 and WriteRays:
    if os.path.exists('ray_files'):
        shutil.rmtree('ray_files')
        os.mkdir('ray_files')
    else:
        os.mkdir('ray_files')


# makes the projections
if comm.rank == 0 and MakeProjections == True:

    for this_scale in proj_scales:

        sp = ds.sphere(c, (this_scale/2.,'kpc'))
    
        for this_dict in field_dicts:
            
            for this_dim in proj_dims:

                proj = yt.ProjectionPlot(ds, 'x', this_dict['field'], weight_field='cell_mass', data_source=sp, width=(this_scale,'kpc'))
                proj.annotate_title(str(ds)+' '+this_dict['title'])
                proj.set_cmap(this_dict['field'], this_dict['colormap'])
                scale = '{:d}'.format(int(this_scale))
                this_filename = str(ds)+'_' + this_dict['file_name'] + '_' + this_dim + '_' + scale + 'kpc.png'
                proj.save(this_filename)
                del proj
        del sp

if Debug:
    comm.Barrier()

# Specify fields to store from the ray
field_list = ['density', 'temperature', 'O_p5_number_density', 'H_p0_number_density', 'metallicity']

# Initialize arrays to use within the loop
H_I_column = np.array([])
O_VI_column = np.array([])


# x is theta, y is phi
Dx = 2.0*np.pi/pixels_per_dim
Dy = np.pi/pixels_per_dim

# Create grid of x & y (theta & phi) values 
x, y = np.mgrid[slice(-np.pi+Dx/2, np.pi+Dx/2, Dx),
                slice(-np.pi/2+Dy/2, np.pi/2+Dy/2, Dy) ]

original_shape = x.shape

if Debug:
    print("Dx, Dy: ", Dx, Dy)
    print("shape is: ", original_shape)

# reshape arrays to use in loop
# so we can divide up rays among processors
x = np.reshape(x,-1)
y = np.reshape(y,-1)

H_I = np.zeros_like(x)
O_VI = np.zeros_like(x)
vel_los = np.zeros_like(x)

# Makes sure rays divide evenly among processors
if (x.size % comm.size == 0):
    if Debug:
        print("Task", comm.rank, "no leftovers!")
else:
    print("Task", comm.rank, "your array size and MPI tasks do not divide evenly!")
    comm.Abort(errorcode=123)

# Dividing the number of rays among the processors
dN = x.size // comm.size
start_index = comm.rank*dN
end_index = (comm.rank+1)*dN

if Debug:
    print("I am on rank", comm.rank, "and my start and end indices are:", start_index, end_index, "of", x.size)

# start the timer for ray generation
raygen_start =  time.time()

for i in range(start_index,end_index):
    theta = x[i]
    phi = y[i]
    
    # convert theta & phi to cartesian coordinates
    dx = R*np.cos(theta)*np.cos(phi)
    dy = R*np.sin(theta)*np.cos(phi)
    dz = R*np.sin(phi)
    
    dv = YTArray([dx, dy, dz], 'kpc')
    
    ray_end = ray_start + dv
    
    padded_num = '{:06d}'.format(i)
    
    if WriteRays:
        rayfile = 'ray_files/ray'+str(ds)+'_'+padded_num+'.h5'
    else:
        rayfile = None
        
    # Create light ray
    lr = trident.LightRay(ds)
    ray = lr.make_light_ray(start_position=ray_start, end_position=ray_end, data_filename=rayfile, fields=field_list)
    
    # Extract data from ray
    ray_data = ray.all_data()
    
    #convert path length to kpc
    path_length = ray_data['dl'].convert_to_units('kpc').d
        
    
    
    # Remove the first N kpc from the ray data
    path = np.zeros(len(path_length))
    
    for h in range(len(path_length)-1):
        dl = path_length[h]
        p = path[h]
        path[h+1] = dl + p
        
    for g in range(len(path)):
        if path[g] < remove_first_N_kpc:
            continue
        else:
            start = g
            break
    
    path_mod = path[start:]
    
    
    # Convert from number density to column density
    H_I_number_density_mod = ray_data['H_p0_number_density'][start:]
    H_I_density = ray_data['dl']*ray_data['H_p0_number_density']
    H_I_density = H_I_density.in_units('cm**-2')
    H_I_density_mod = H_I_density[start:]
    H_I_column_density = sum(H_I_density_mod.d)
    H_I_column = np.append(H_I_column, H_I_column_density)        
        
    O_VI_number_density_mod = ray_data['O_p5_number_density'][start:]
    O_VI_density = ray_data['dl']*ray_data['O_p5_number_density']
    O_VI_density = O_VI_density.in_units('cm**-2')
    O_VI_density_mod = O_VI_density[start:]
    O_VI_column_density = sum(O_VI_density_mod.d)
    O_VI_column = np.append(O_VI_column, O_VI_column_density)
    
    # Convert LOS velocity to km/s weighted by the HI density and calculate the mean
    los_mod = ray_data['velocity_los'][start:]
    los_mod = los_mod.in_units('km/s')
    
    los_mod = H_I_density_mod*los_mod
    mean_los = sum(los_mod.d)/sum(H_I_density_mod.d)
    
    
    # Save column densities & mean velocity to arrays
    H_I[i] += H_I_column_density
    O_VI[i] += O_VI_column_density
    vel_los[i] += mean_los
    
    # Make plots of various quantities along the line of sight
    r = npr.uniform()
    if r <= fraction_ray_plots:
        plt.subplot(321)
        plt.semilogy(path_mod, H_I_number_density_mod, color='purple')
        #plt.ylim(10**-12.5, 10**-9.5)
        plt.xlabel("Distance (kpc)")
        plt.ylabel('H I Density')
        plt.title("H I density")
        plt.grid()        

        plt.subplot(322)
        plt.semilogy(path_mod, O_VI_number_density_mod, color='blue')
        #plt.ylim(10**-16, 10**-8)
        plt.xlabel("Distance (kpc)")
        plt.ylabel('O VI Density')
        plt.title("O VI density")
        plt.grid()
       
        plt.subplot(323)
        plt.semilogy(path, ray_data['temperature'], 'r-')
        #plt.ylim(10**5, 10**8.5)
        plt.title('Temperature')
        plt.xlabel('Distance (kpc)')
        plt.ylabel('Temperature (K)')
        plt.grid()
    
        plt.subplot(324)
        plt.semilogy(path, ray_data['density'], 'g-')
        #plt.ylim(10**-28.8, 10**-26)
        plt.xlabel("Distance (kpc)")
        plt.ylabel('Density (g/cm^2)')
        plt.title("Density")
        plt.grid()
	
        plt.subplot(325)
        plt.semilogy(path, ray_data['metallicity'], linestyle='-',color='orange')
        plt.xlabel('Distance (kpc)')
        plt.ylabel('Metallicity')
        plt.title('Metallicity')
        plt.grid()

        v_los = ray_data['velocity_los'].in_units('km/s')        

        plt.subplot(326)
        plt.plot(path, v_los, color='black')
        plt.xlabel('Distance (kpc)')
        plt.ylabel('Velocity (km/s)')
        plt.title('Line of Sight Velocity')
        plt.grid()

        plt.tight_layout()
        plt.savefig('ray_data_'+str(ds)+'_'+str(i)+'.png')
       
        plt.clf()

# end timer for ray generation
raygen_end =  time.time()

time_per_ray = np.array([(raygen_end-raygen_start)/(end_index-start_index)])

tpr_min = np.array([-1.0])
tpr_max = np.array([-1.0])
tpr_mean = np.array([-1.0])
#tpr_min = tmr_max = tpr_mean = -1.0

# get some communication information
comm.Reduce(time_per_ray,tpr_min,op=MPI.MIN)
comm.Reduce(time_per_ray,tpr_max,op=MPI.MAX)
comm.Reduce(time_per_ray,tpr_mean,op=MPI.SUM)
tpr_mean /= comm.size

reduce_HI = np.zeros_like(H_I)
reduce_OVI = np.zeros_like(O_VI)
reduce_vlos = np.zeros_like(vel_los)

# Combines arrays from each processor
comm.Reduce(H_I,reduce_HI,op=MPI.SUM)
comm.Reduce(O_VI,reduce_OVI,op=MPI.SUM)
comm.Reduce(vel_los,reduce_vlos,op=MPI.SUM)



# Save arrays into a text file on the main processor
if comm.rank == 0:
    filename = str(ds)+'_data.txt'
    np.savetxt(filename,np.c_[x, y, reduce_HI, reduce_OVI, reduce_vlos],fmt="%5.4e",header="theta, phi, HI, OVI, V_LOS")


x = np.reshape(x,original_shape)
y = np.reshape(y,original_shape)
H_I = np.reshape(reduce_HI,original_shape)
O_VI = np.reshape(reduce_OVI,original_shape)
vel_los = np.reshape(reduce_vlos,original_shape)

if Debug:
    print("Ray generation took {:.3f} seconds on MPI task {:d}".format(raygen_end-raygen_start, comm.rank))
    print("That's {:.3e} seconds per ray on this MPI task".format(time_per_ray[0]))


if comm.rank == 0:
    plt.figure()

    # sets up the plot as an Aitoff projection
    plt.subplot(111, projection="aitoff")

    # Convert to log scale
    logH_I = np.log10(H_I + 10.0**-20.0)

    # Aitoff projection using pcolor
    plt.pcolor(x, y, logH_I, cmap='plasma') #, vmin=12.5, vmax=21.5)
    cbar=plt.colorbar(pad=0.02,shrink=0.55)
    cbar.ax.tick_params(labelsize=6)
    cbar.set_label(r'log N$_{HI}$ [cm$^{-2}$]',size=8)
    plt.title(str(ds)+" HI Column Density",y=1.05,size=10)
    plt.grid(True,alpha=0.5)
    plt.xticks(fontsize=7)
    plt.yticks(fontsize=7)
    plt.savefig('HI_aitoff_'+str(ds)+'.png',bbox_inches='tight',dpi=400)

    plt.clf()


    # sets up the plot as an Aitoff projection
    plt.subplot(111, projection="aitoff")

    # Convert to log scale
    logO_VI = np.log10(O_VI + 10**-20.0)

    # Aitoff projection using pcolor
    plt.pcolor(x, y, logO_VI, cmap='viridis') #, vmin=12.0, vmax=16.0)
    plt.title(str(ds)+" OVI Column Density",y=1.05,size=10)
    plt.grid(True,alpha=0.5)
    cbar = plt.colorbar(pad=0.02,shrink=0.55)
    cbar.ax.tick_params(labelsize=6)
    cbar.set_label(r'N$_{OVI}$ [cm$^{-2}$]',size=8)
    plt.xticks(fontsize=7)
    plt.yticks(fontsize=7)
    plt.savefig('OVI_aitoff_'+str(ds)+'.png',bbox_inches='tight',dpi=400)

    plt.clf()


    # sets up the plot as an Aitoff projection
    plt.subplot(111, projection="aitoff")

    # Aitoff projection using pcolor
    plt.pcolor(x, y, vel_los, cmap='plasma') #, vmin=-300.0, vmax=+300.0)
    cbar=plt.colorbar(pad=0.02,shrink=0.55)
    cbar.ax.tick_params(labelsize=6)
    cbar.set_label(r'v$_{LOS}$ [km/s]',size=8)
    plt.title(str(ds)+" HI line-of-sight velocity",y=1.05,size=10)
    plt.grid(True,alpha=0.5)
    plt.xticks(fontsize=7)
    plt.yticks(fontsize=7)
    plt.savefig('HI_LOS_aitoff_'+str(ds)+'.png',bbox_inches='tight',dpi=400)

    
    plt.clf()

    # Repeat projections, but with Cartesian projection instead of Aitoff
    plt.subplot(111)

    plt.pcolor(x, y, logO_VI, cmap='viridis')
    plt.title(str(ds)+" OVI Column Density",y=1.05,size=10)
    plt.grid(True,alpha=0.5)
    cbar = plt.colorbar(pad=0.02,shrink=0.55)
    cbar.ax.tick_params(labelsize=6)
    cbar.set_label(r'N$_{OVI}$ [cm$^{-2}$]',size=8)
    plt.xticks(fontsize=7)
    plt.yticks(fontsize=7)
    plt.axis([x.min(), x.max(), y.min(), y.max()])
    plt.savefig('OVI_cartesian_'+str(ds)+'.png',bbox_inches='tight',dpi=400)

    plt.clf()

    plt.subplot(111)

    plt.pcolor(x, y, logH_I, cmap='plasma')
    cbar=plt.colorbar(pad=0.02,shrink=0.55)
    cbar.ax.tick_params(labelsize=6)
    cbar.set_label(r'log N$_{HI}$ [cm$^{-2}$]',size=8)
    plt.title(str(ds)+" HI Column Density",y=1.05,size=10)
    plt.grid(True,alpha=0.5)
    plt.xticks(fontsize=7)
    plt.yticks(fontsize=7)
    plt.axis([x.min(), x.max(), y.min(), y.max()])
    plt.savefig('HI_cartesian_'+str(ds)+'.png',bbox_inches='tight',dpi=400)

    plt.clf()

    plt.pcolor(x, y, vel_los, cmap='plasma')
    cbar=plt.colorbar(pad=0.02,shrink=0.55)
    cbar.ax.tick_params(labelsize=6)
    cbar.set_label(r'v$_{LOS}$ [km/s]',size=8)
    plt.title(str(ds)+" HI line-of-sight velocity",y=1.05,size=10)
    plt.grid(True,alpha=0.5)
    plt.xticks(fontsize=7)
    plt.yticks(fontsize=7)
    plt.axis([x.min(), x.max(), y.min(), y.max()])
    plt.savefig('HI_LOS_cartesian_'+str(ds)+'.png',bbox_inches='tight',dpi=400)

    plt.clf()

    logH_I = np.reshape(logH_I, -1)
    logO_VI = np.reshape(logO_VI, -1)

    # Create plot of cumulative covering fraction for H I and O VI
    plt.subplot(111)
    logH_I = np.sort(logH_I)
    f = np.array([])
    for a in range(len(logH_I)):
        num = logH_I[a:].size/logH_I.size
        f = np.append(f, num)


    logf = np.log10(f + 10**-20)

    plt.plot(logH_I, logf, 'bo-')
    plt.xlabel('log HI')
    plt.ylabel('log f(>N)')
    plt.title('HI cumulative covering fraction')
    plt.savefig(str(ds)+'_HI_column_dens_fraction.png',bbox_inches='tight',dpi=400)

    plt.clf()


    plt.subplot(111)
    logO_VI = np.sort(logO_VI)
    g = np.array([])
    for a in range(len(logO_VI)):
        num = logO_VI[a:].size/logO_VI.size
        g = np.append(g, num)


    logg = np.log10(g + 10**-20)

    plt.plot(logO_VI, logg, 'ro-')
    plt.xlabel('log OVI')
    plt.ylabel('log f(>N)')
    plt.title('OVI cumulative covering fraction')
    plt.savefig(str(ds)+'_OVI_column_dens_fraction.png',bbox_inches='tight',dpi=400)


comm.Barrier()

t_end = time.time()

if comm.rank == 0:
    print("This calculation took {:.3f} seconds".format(t_end-t_start))
    print("That's {:.3e} seconds per ray per MPI task (including serial operations)".format( (t_end-t_start)*comm.size/reduce_HI.size))
    print("For just ray generation, mean (min/max) time per MPI task is:  {:.3e} ({:.3e}/{:.3e})".format(tpr_mean[0],tpr_min[0],tpr_max[0]))


if comm.rank == 0:
    print("JOB FINISHED")
