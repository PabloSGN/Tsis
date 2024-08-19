import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

def plot_polar_azimuthal(lat_range, lon_range, projection='npstere', bounding_lat=60):
    """
    Plot a polar map with a given latitude and longitude range using a Polar Azimuthal Equidistant Projection.
    
    Parameters:
    lat_range (tuple): The range of latitudes to display (min_lat, max_lat).
    lon_range (tuple): The range of longitudes to display (min_lon, max_lon).
    projection (str): The type of polar projection ('npstere' for North Pole, 'spstere' for South Pole).
    bounding_lat (float): The bounding latitude for the polar stereographic projection.
    """
    plt.figure(figsize=(10, 10))
    
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
    
    plt.title(f'Polar Azimuthal Equidistant Projection\nLatitude range: {lat_range}, Longitude range: {lon_range}')

# Example: Plot North Polar region with latitude range 60 to 90 and longitude range -180 to 180
plot_polar_azimuthal(lat_range=(10, 90), lon_range=(-180, 180), projection='npstere', bounding_lat=10)



plt.show()



