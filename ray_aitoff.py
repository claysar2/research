import matplotlib as mpl
mpl.use('Agg')

import yt
import trident
import random
import numpy as np
from yt import YTArray
import matplotlib.pyplot as plt
import math
import numpy.random as npr


# Load dataset
fn = "/mnt/scratch/dsilvia/simulations/reu_sims/MW_1638kpcBox_800pcCGM_200pcDisk_lowres/DD1500/DD1500"
ds = yt.load(fn)

# Add H I & O VI ion fields using Trident
trident.add_ion_fields(ds, ions=['O VI', 'H I'], ftype="gas")

# Specify position where the ray starts
c = ds.arr([0.5, 0.5, 0.5], 'unitary')
c = c.in_units('kpc')
X = YTArray([8., 0., 0.], 'kpc')
ray_start = c - X

# Length of Ray
R = 200.0


# Create data spheres of radii 15 kpc and 150 kpc
sp1 = ds.sphere(c, (15.,'kpc'))
sp2 = ds.sphere(c, (150.,'kpc'))

dens1x = yt.ProjectionPlot(ds, 'x', 'density', weight_field='cell_mass', data_source=sp1, width=(50.,'kpc'))
dens1x.annotate_title(str(ds)+' Density')
dens1x.set_cmap('density', 'kelp')
dens1x.save(str(ds)+'_density_x_15kpc.png')
dens1y = yt.ProjectionPlot(ds, 'y', 'density', weight_field='cell_mass', data_source=sp1, width=(50.,'kpc'))
dens1y.annotate_title(str(ds)+' Density')
dens1y.set_cmap('density', 'kelp')
dens1y.save(str(ds)+'_density_y_15kpc.png')
dens1z = yt.ProjectionPlot(ds, 'z', 'density', weight_field='cell_mass', data_source=sp1, width=(50.,'kpc'))
dens1z.annotate_title(str(ds)+' Density')
dens1z.set_cmap('density', 'kelp')
dens1z.save(str(ds)+'_density_z_15kpc.png')
dens2x = yt.ProjectionPlot(ds, 'x', 'density', weight_field='cell_mass', data_source=sp2, width=(300.,'kpc'))
dens2x.annotate_title(str(ds)+' Density')
dens2x.set_cmap('density', 'kelp')
dens2x.save(str(ds)+'_density_x_150kpc.png')
dens2y = yt.ProjectionPlot(ds, 'y', 'density', weight_field='cell_mass', data_source=sp2, width=(300.,'kpc'))
dens2y.annotate_title(str(ds)+' Density')
dens2y.set_cmap('density', 'kelp')
dens2y.save(str(ds)+'_density_y_150kpc.png')
dens2z = yt.ProjectionPlot(ds, 'z', 'density', weight_field='cell_mass', data_source=sp2, width=(300.,'kpc'))
dens2z.annotate_title(str(ds)+' Density')
dens2z.set_cmap('density', 'kelp')
dens2z.save(str(ds)+'_density_z_150kpc.png')

temp1x = yt.ProjectionPlot(ds, 'x', 'temperature', weight_field='cell_mass', data_source=sp1, width=(50.,'kpc'))
temp1x.annotate_title(str(ds)+' Temperature')
temp1x.set_cmap('temperature', 'RdBu_r')
temp1x.save(str(ds)+'_temp_x_15kpc.png')
temp1y = yt.ProjectionPlot(ds, 'y', 'temperature', weight_field='cell_mass', data_source=sp1, width=(50.,'kpc'))
temp1y.annotate_title(str(ds)+' Temperature')
temp1y.set_cmap('temperature', 'RdBu_r')
temp1y.save(str(ds)+'_temp_y_15kpc.png')
temp1z = yt.ProjectionPlot(ds, 'z', 'temperature', weight_field='cell_mass', data_source=sp1, width=(50.,'kpc'))
temp1z.annotate_title(str(ds)+' Temperature')
temp1z.set_cmap('temperature', 'RdBu_r')
temp1z.save(str(ds)+'_temp_z_15kpc.png')
temp2x = yt.ProjectionPlot(ds, 'x', 'temperature', weight_field='cell_mass', data_source=sp2, width=(300.,'kpc'))
temp2x.annotate_title(str(ds)+' Temperature')
temp2x.set_cmap('temperature', 'RdBu_r')
temp2x.save(str(ds)+'_temp_x_150kpc.png')
temp2y = yt.ProjectionPlot(ds, 'y', 'temperature', weight_field='cell_mass', data_source=sp2, width=(300.,'kpc'))
temp2y.annotate_title(str(ds)+' Temperature')
temp2y.set_cmap('temperature', 'RdBu_r')
temp2y.save(str(ds)+'_temp_y_150kpc.png')
temp2z = yt.ProjectionPlot(ds, 'z', 'temperature', weight_field='cell_mass', data_source=sp2, width=(300.,'kpc'))
temp2z.annotate_title(str(ds)+' Temperature')
temp2z.set_cmap('temperature', 'RdBu_r')
temp2z.save(str(ds)+'_temp_z_150kpc.png')

meta1x = yt.ProjectionPlot(ds, 'x', 'metallicity', weight_field='cell_mass', data_source=sp1, width=(50.,'kpc'))
meta1x.annotate_title(str(ds)+' Metallicity')
meta1x.set_cmap('metallicity', 'arbre')
meta1x.save(str(ds)+'_metallicity_x_15kpc.png')
meta1y = yt.ProjectionPlot(ds, 'y', 'metallicity', weight_field='cell_mass', data_source=sp1, width=(50.,'kpc'))
meta1y.annotate_title(str(ds)+' Metallicity')
meta1y.set_cmap('metallicity', 'arbre')
meta1y.save(str(ds)+'_metallicity_y_15kpc.png')
meta1z = yt.ProjectionPlot(ds, 'z', 'metallicity', weight_field='cell_mass', data_source=sp1, width=(50.,'kpc'))
meta1z.annotate_title(str(ds)+' Metallicity')
meta1z.set_cmap('metallicity', 'arbre')
meta1z.save(str(ds)+'_metallicity_z_15kpc.png')
meta2x = yt.ProjectionPlot(ds, 'x', 'metallicity', weight_field='cell_mass', data_source=sp2, width=(300.,'kpc'))
meta2x.annotate_title(str(ds)+' Metallicity')
meta2x.set_cmap('metallicity', 'arbre')
meta2x.save(str(ds)+'_metallicity_x_150kpc.png')
meta2y = yt.ProjectionPlot(ds, 'y', 'metallicity', weight_field='cell_mass', data_source=sp2, width=(300.,'kpc'))
meta2y.annotate_title(str(ds)+' Metallicity')
meta2y.set_cmap('metallicity', 'arbre')
meta2y.save(str(ds)+'_metallicity_y_150kpc.png')
meta2z = yt.ProjectionPlot(ds, 'z', 'metallicity', weight_field='cell_mass', data_source=sp2, width=(300.,'kpc'))
meta2z.annotate_title(str(ds)+' Metallicity')
meta2z.set_cmap('metallicity', 'arbre') 
meta2z.save(str(ds)+'_metallicity_z_150kpc.png')

densHI1x = yt.ProjectionPlot(ds, 'x', 'H_p0_number_density', weight_field='cell_mass', data_source=sp1, width=(50.,'kpc'))
densHI1x.annotate_title(str(ds)+' H I Density')
densHI1x.set_cmap('H_p0_number_density', 'plasma')
densHI1x.save(str(ds)+'_H_I_density_x_15kpc.png')
densHI1y = yt.ProjectionPlot(ds, 'y', 'H_p0_number_density', weight_field='cell_mass', data_source=sp1, width=(50.,'kpc'))
densHI1y.annotate_title(str(ds)+' H I Density')
densHI1y.set_cmap('H_p0_number_density', 'plasma')
densHI1y.save(str(ds)+'_H_I_density_y_15kpc.png')
densHI1z = yt.ProjectionPlot(ds, 'z', 'H_p0_number_density', weight_field='cell_mass', data_source=sp1, width=(50.,'kpc'))
densHI1z.annotate_title(str(ds)+' H I Density')
densHI1z.set_cmap('H_p0_number_density', 'plasma')
densHI1z.save(str(ds)+'_H_I_density_z_15kpc.png')
densHI2x = yt.ProjectionPlot(ds, 'x', 'H_p0_number_density', weight_field='cell_mass', data_source=sp2, width=(300.,'kpc'))
densHI2x.annotate_title(str(ds)+' H I Density')
densHI2x.set_cmap('H_p0_number_density', 'plasma')
densHI2x.save(str(ds)+'_H_I_density_x_150kpc.png')
densHI2y = yt.ProjectionPlot(ds, 'y', 'H_p0_number_density', weight_field='cell_mass', data_source=sp2, width=(300.,'kpc'))
densHI2y.annotate_title(str(ds)+' H I Density')
densHI2y.set_cmap('H_p0_number_density', 'plasma')
densHI2y.save(str(ds)+'_H_I_density_y_150kpc.png')
densHI2z = yt.ProjectionPlot(ds, 'z', 'H_p0_number_density', weight_field='cell_mass', data_source=sp2, width=(300.,'kpc'))
densHI2z.annotate_title(str(ds)+' H I Density')
densHI2z.set_cmap('H_p0_number_density', 'plasma')
densHI2z.save(str(ds)+'_H_I_density_z_150kpc.png')

densOVI1x = yt.ProjectionPlot(ds, 'x', 'O_p5_number_density', weight_field='cell_mass', data_source=sp1, width=(50.,'kpc'))
densOVI1x.annotate_title(str(ds)+' O VI Density')
densOVI1x.set_cmap('O_p5_number_density', 'viridis')
densOVI1x.save(str(ds)+'_O_VI_density_x_15kpc.png')
densOVI1y = yt.ProjectionPlot(ds, 'y', 'O_p5_number_density', weight_field='cell_mass', data_source=sp1, width=(50.,'kpc'))
densOVI1y.annotate_title(str(ds)+' O VI Density')
densOVI1y.set_cmap('O_p5_number_density', 'viridis')
densOVI1y.save(str(ds)+'_O_VI_density_y_15kpc.png')
densOVI1z = yt.ProjectionPlot(ds, 'z', 'O_p5_number_density', weight_field='cell_mass', data_source=sp1, width=(50.,'kpc'))
densOVI1z.annotate_title(str(ds)+' O VI Density')
densOVI1z.set_cmap('O_p5_number_density', 'viridis')
densOVI1z.save(str(ds)+'_O_VI_density_z_15kpc.png')
densOVI2x = yt.ProjectionPlot(ds, 'x', 'O_p5_number_density', weight_field='cell_mass', data_source=sp2, width=(300.,'kpc'))
densOVI2x.annotate_title(str(ds)+' O VI Density')
densOVI2x.set_cmap('O_p5_number_density', 'viridis')
densOVI2x.save(str(ds)+'_O_VI_density_x_150kpc.png')
densOVI2y = yt.ProjectionPlot(ds, 'y', 'O_p5_number_density', weight_field='cell_mass', data_source=sp2, width=(300.,'kpc'))
densOVI2y.annotate_title(str(ds)+' O VI Density')
densOVI2y.set_cmap('O_p5_number_density', 'viridis')
densOVI2y.save(str(ds)+'_O_VI_density_y_150kpc.png')
densOVI2z = yt.ProjectionPlot(ds, 'z', 'O_p5_number_density', weight_field='cell_mass', data_source=sp2, width=(300.,'kpc'))
densOVI2z.annotate_title(str(ds)+' O VI Density')
densOVI2z.set_cmap('O_p5_number_density', 'viridis')
densOVI2z.save(str(ds)+'_O_VI_density_z_150kpc.png')


# Specify fields to store from the ray
field_list = ['density', 'temperature', 'O_p5_number_density', 'H_p0_number_density', 'metallicity']

# Initialize arrays to use within the loop
H_I_column = np.array([])
O_VI_column = np.array([])

# x is theta, y is phi
Dx, Dy = 0.1, 0.05
x, y = np.mgrid[slice(-np.pi, np.pi + Dx, Dx),
                slice(-np.pi/2, np.pi/2 + Dy, Dy) ]
# Create H I & O VI column density arrays 
H_I = np.zeros(x.shape)
O_VI = np.zeros(x.shape)
vel_los = np.zeros(x.shape)

a = np.arange(x.shape[0])
b = np.arange(x.shape[1])

for i in a:
    for j in b:
        
        theta = x[i,j]
        phi = y[i,j]
    
        dx = R*np.cos(theta)*np.sin(phi)
        dy = R*np.sin(theta)*np.sin(phi)
        dz = R*np.cos(phi)

        dv = YTArray([dx, dy, dz], 'kpc')
    
        ray_end = ray_start + dv
       
        rayfile = 'ray_files/ray'+str(ds)+'_'+str(i)+'_'+str(j)+'.h5'
    
        ray = trident.make_simple_ray(ds,
                                      start_position=ray_start,
                                      end_position=ray_end,
                                      data_filename=rayfile,
                                      fields=field_list,
                                      ftype='gas')
    
        ray_data = ray.all_data()

        path_length = ray_data['dl'].convert_to_units('kpc').d
        


       # Remove the first 5 kpc from the ray data
        path = np.zeros(len(path_length))
        
        for h in range(len(path_length)-1):
            dl = path_length[h]
            p = path[h]
            path[h+1] = dl + p
        
        for g in range(len(path)):
            if path[g] < 5:
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
        

        # Save column densities to arrays
        H_I[i,j] += H_I_column_density
        O_VI[i,j] += O_VI_column_density
	
        r = np.random.uniform()

        if r <= 0.01:
            plt.subplot(321)
            plt.semilogy(path_mod, H_I_number_density_mod, color='purple')
            plt.ylim(10**-12.5, 10**-9.5)
            plt.xlabel("Distance (kpc)")
            plt.ylabel('H I Density')
            plt.title("H I density")
            plt.grid()        

            plt.subplot(322)
            plt.semilogy(path_mod, O_VI_number_density_mod, color='blue')
            plt.ylim(10**-16, 10**-8)
            plt.xlabel("Distance (kpc)")
            plt.ylabel('O VI Density')
            plt.title("O VI density")
            plt.grid()
       
            plt.subplot(323)
            plt.semilogy(path, ray_data['temperature'], 'r-')
            plt.ylim(10**5, 10**8.5)
            plt.title('Temperature')
            plt.xlabel('Distance (kpc)')
            plt.ylabel('Temperature (K)')
            plt.grid()
    
            plt.subplot(324)
            plt.semilogy(path, ray_data['density'], 'g-')
            plt.ylim(10**-28.8, 10**-26)
            plt.xlabel("Distance (kpc)")
            plt.ylabel('Density (g/cm^2)')
            plt.title("Density")
            plt.grid()
	
            plt.subplot(325)
            plt.semilogy(path, ray_data['metallicity'], color='orange')
            plt.xlabel('Distance (kpc)')
            plt.ylabel('Metallicity')
            plt.title('Metallicity')
            plt.grid()

            v_los = ray_data['velocity_los'].in_units('km/s')
        
            plt.subplot(326)
            plt.semilogy(path, v_los, color='black')
            plt.xlabel('Distance (kpc)')
            plt.ylabel('Velocity (km/s)')
            plt.title('Line of Sight Velocity')
            plt.grid()

            plt.tight_layout()
            plt.savefig('ray_data_'+str(ds)+'_'+str(i)+'_'+str(j)+'.png')


plt.figure()

# sets up the plot as an Aitoff projection
plt.subplot(111, projection="aitoff")

# Convert to log scale
logH_I = np.log10(H_I + 10.0**-20.0)

# Aitoff projection using pcolor
plt.pcolor(x, y, logH_I, cmap='plasma', vmin=12.5, vmax=21.5)
plt.colorbar(label='cm^-2')
plt.title(str(ds)+" H I Column Density")
plt.grid(True)
plt.savefig('HI_aitoff_'+str(ds)+'.png')


plt.clf()


# sets up the plot as an Aitoff projection
plt.subplot(111, projection="aitoff")

# Convert to log scale
logO_VI = np.log10(O_VI + 10**-20.0)

# Aitoff projection using pcolor
plt.pcolor(x, y, logO_VI, cmap='viridis', vmin=12.0, vmax=16.0)
plt.title(str(ds)+" O VI Column Density")
plt.grid(True)
plt.colorbar(label='cm^-2')
plt.savefig('OVI_aitoff_'+str(ds)+'.png')

plt.clf()


plt.subplot(111)
logH_I = np.sort(logH_I)
f = np.array([])
for a in range(len(logH_I)):
    num = logH_I[a:].size/logH_I.size
    f = np.append(f, num)


logf = np.log10(f + 10**-20)

plt.plot(logH_I, logf, 'o')
plt.xlabel('log H I')
plt.ylabel('log N')
plt.savefig(str(ds)+'_HI_column_dens.png')

plt.clf()


plt.subplot(111)
logO_VI = np.sort(logO_VI)
g = np.array([])
for a in range(len(logO_VI)):
    num = logO_VI[a:].size/logO_VI.size
    g = np.append(g, num)


logg = np.log10(g + 10**-20)

plt.plot(logO_VI, logg, 'o')
plt.xlabel('log	O VI')
plt.ylabel('log	N')
plt.savefig(str(ds)+'_OVI_column_dens.png')



# reshape data to 1d arrays for saving into txt file
H_I = np.reshape(H_I, -1)
O_VI = np.reshape(O_VI, -1)
x = np.reshape(x, -1)
y = np.reshape(y, -1)


# Save array data into text file
filename = str(ds)+'_data.txt'
np.savetxt(filename,np.c_[x,y,H_I,O_VI],fmt="%5.4f",header="theta, phi, H I, O VI")

