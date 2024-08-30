import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import matplotlib as mp
import numpy as np


plt.rc('font', family='serif')
plt.rc('xtick', labelsize='small')
plt.rc('ytick', labelsize='small')
mp.rcParams["font.size"] = "12.5"

def plot_polar_azimuthal(lat_range, lon_range, trajectory, projection='npstere', bounding_lat=60):
    """
    Plot a polar map with a given latitude and longitude range using a Polar Azimuthal Equidistant Projection.
    
    Parameters:
    lat_range (tuple): The range of latitudes to display (min_lat, max_lat).
    lon_range (tuple): The range of longitudes to display (min_lon, max_lon).
    projection (str): The type of polar projection ('npstere' for North Pole, 'spstere' for South Pole).
    bounding_lat (float): The bounding latitude for the polar stereographic projection.
    """
    plt.figure(figsize=(5.5, 5.5))
    
    # Set up the Basemap
    if projection == 'npstere':
        m = Basemap(projection=projection, boundinglat=bounding_lat, lon_0=0, resolution='l')
    elif projection == 'spstere':
        m = Basemap(projection=projection, boundinglat=-bounding_lat, lon_0=0, resolution='l')
    else:
        raise ValueError("Projection must be 'npstere' for North Pole or 'spstere' for South Pole.")
    
    m.drawcoastlines()
    m.drawcountries()
    m.fillcontinents(color='mediumseagreen', lake_color='turquoise')
    m.drawmapboundary(fill_color='turquoise')
    m.drawparallels(range(int(lat_range[0]), int(lat_range[1]), 10))
    m.drawmeridians(range(int(lon_range[0]), int(lon_range[1]), 30))
    
    x, y = m(trajectory[:, 1], trajectory[:, 0])

    m.plot(x, y, color = 'crimson', lw = 4)

    plt.title("Sunrise III  trajectory.")

    return m

# Example: Plot North Polar region with latitude range 60 to 90 and longitude range -180 to 180

coords = np.loadtxt("Track740N_mod.txt", usecols = (0, 1))

print(np.shape(coords))
print(coords[:, 1])

mapp = plot_polar_azimuthal(lat_range=(10, 90), lon_range=(-180, 180), trajectory=coords, projection='npstere', bounding_lat=10)


labels = [r"Esrange, $10^{th}$ of June 2024", r"Landing Point, $16^{th}$ of June 2024"]

lab_lat = [coords[0, 0], coords[-1, 0]]
lab_lon = [coords[0, 1], coords[-1, 1]]

colors = ["indigo", "darkorange"]

for i, label in enumerate(labels):

    x, y = mapp(lab_lon[i], lab_lat[i])
    
    mapp.plot(x, y, color = colors[i], marker = 'X', markersize = 10, ls = '', label = labels[i])
    


plt.legend(edgecolor = 'k')




plt.show()


