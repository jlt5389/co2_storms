#!/usr/bin/python
# Storm Systems Effect On CO2
#   Courtney McDaniels
#   Tyanna Hellams
#   Joshua Tscheschlog
# Version: 07/08/2016

# === Import =================
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from datetime import datetime, timedelta
from mpl_toolkits.basemap import Basemap, addcyclic
from scipy.ndimage.filters import minimum_filter, maximum_filter
from netCDF4 import Dataset
import math as ma
import webbrowser

# === Global Variables  ======
start_t = 740   # Start Date
end_t = 780     # End Date
storm_sensitivity = 12  # Lower is more sensitive
filename = '/home/meteo/cmm6461/meteo473/spring2015/milestone4/mslp.2010.nc'
data = Dataset(filename)  # Open OpenDAP dataset.
e_rad = 6371    # Radius of Earth (km)
mini_dist = 500 # How close a storm must be in 2 different periods to be considered the same storm (km)
lats = data.variables['lat'][:] # Read Latitudes From File
lons = data.variables['lon'][:] # Read Longitudes From File
strm_case = '93'   # Choose a storm case study to study changes in CO2 around storm area
co2_d = 10       # How many days before/after to mean co2 data
co2_area = 3    # How many degrees around the storm path to average co2 levels for analysis

# === Functions ==============
def extrema(mat,mode='wrap',window=5):
    """find the indices of local extrema (min and max)
    in the input array."""
    mn = minimum_filter(mat, size=window, mode=mode)
    mx = maximum_filter(mat, size=window, mode=mode)
    # (mat == mx) true if pixel is equal to the local max
    # (mat == mn) true if pixel is equal to the local in
    # Return the indices of the maxima, minima
    return np.nonzero(mat == mn), np.nonzero(mat == mx)

# Function that converts the native date/time index to a string in format "yyyy-mm-dd hh UTC"
def hour_date(time_list, hour_index):
    d = datetime(1800, 1, 1)
    d = d + timedelta(hours = time_list[hour_index])
    return d

# Degrees to Radians
def to_rad(degs):
    radss = degs*3.1415/180
    return radss

# Returns distance between two coordinates on Earth's surface (km)
def great_circle(lat_t1, lat_t2, lon_t1, lon_t2):
    lat_t1 = to_rad(lat_t1)
    lat_t2 = to_rad(lat_t2)
    lon_t1 = to_rad(lon_t1)
    lon_t2 = to_rad(lon_t2)
    if (ma.cos(lat_t1)*ma.cos(lat_t2)*ma.cos(lon_t1-lon_t2)+ma.sin(lat_t1)*ma.sin(lat_t2)) > -1 and (ma.cos(lat_t1)*ma.cos(lat_t2)*ma.cos(lon_t1-lon_t2)+ma.sin(lat_t1)*ma.sin(lat_t2)) < 1:
        e_dist = e_rad*(ma.acos(ma.cos(lat_t1)*ma.cos(lat_t2)*ma.cos(lon_t1-lon_t2)+ma.sin(lat_t1)*ma.sin(lat_t2)))
    else:
        e_dist = 1000
    return e_dist

 
# function to create Basemap instance.
def b_map(lon0, lat0, lon1, lat1, proj):
    b = Basemap(llcrnrlon=0,llcrnrlat=-80,urcrnrlon=360,urcrnrlat=80,projection='mill')
    return b

# STEP 1: Initialize a 'master' data structure and populate it with storm data from the first selected period ==========

# read prmsl for Jan 1, convert to hPa (mb).
# the window parameter controls the number of highs and lows detected.
# (higher value, fewer highs and lows)
master_prmsl = []
master_lmin = []
master_t = []
counter = 0
for date_t in range(start_t,end_t):
    master_prmsl.append(0.01*np.array(data.variables['mslp'][:])[date_t])
    local_min, local_max = extrema(master_prmsl[counter], mode='wrap', window=storm_sensitivity)
    master_lmin.append(local_min)
    master_t.append(data.variables['time'][date_t])
    counter += 1
for date_t in range(end_t+1,end_t+(co2_d*4)):
    master_t.append(data.variables['time'][date_t])
for date_t in range(start_t-(co2_d*4),start_t-1):
    master_t.append(data.variables['time'][date_t])

# Create storm_stats dictionary, populate with all storms at time = 0
storm_stats = {}
for z in range(len(master_lmin[0][0])):
    storm_stats[str(z)] = {}
    storm_stats[str(z)]['lats'] = [lats[master_lmin[0][0][z]]]
    storm_stats[str(z)]['lons'] = [lons[master_lmin[0][1][z]]]
    storm_stats[str(z)]['time'] = [0]
    storm_stats[str(z)]['prmsl'] = [master_prmsl[0][master_lmin[0]][z]]

# STEP 2: Update data structure, organized by individual storm

# Compare each storm for each period with all storms from the previous period in time and measure distance, recording the smallest value.
# If this value exceeds "mini-dist" km, append storm_stats with a new storm, otherwise append closest existing storm
for date_t in range(1,(end_t-start_t)):
    for q in range(len(master_lmin[date_t][0])):
        c_gs = mini_dist # Reset for next storm-check       
        for q2 in range(len(storm_stats)):
            if date_t-1 == storm_stats[str(q2)]['time'][-1]:
                gs = great_circle(lats[master_lmin[date_t][0][q]], storm_stats[str(q2)]['lats'][-1], lons[master_lmin[date_t][1][q]], storm_stats[str(q2)]['lons'][-1])
                if (gs < c_gs):
                    c_gs = gs
                    min_gsi = q2
        # Add new storm
        if c_gs >= mini_dist:
            storm_stats[str(len(storm_stats))] = {}
            storm_stats[str(len(storm_stats)-1)]['lats'] = [lats[master_lmin[date_t][0][q]]]
            storm_stats[str(len(storm_stats)-1)]['lons'] = [lons[master_lmin[date_t][1][q]]]
            storm_stats[str(len(storm_stats)-1)]['time'] = [date_t]
            storm_stats[str(len(storm_stats)-1)]['prmsl'] = [master_prmsl[date_t][master_lmin[date_t]][q]]
        # Update existing storm
        elif c_gs < mini_dist and date_t == storm_stats[str(min_gsi)]['time'][-1] + 1:
            storm_stats[str(min_gsi)]['lats'].append(lats[master_lmin[date_t][0][q]])
            storm_stats[str(min_gsi)]['lons'].append(lons[master_lmin[date_t][1][q]])
            storm_stats[str(min_gsi)]['time'].append(date_t)
            storm_stats[str(min_gsi)]['prmsl'].append(master_prmsl[date_t][master_lmin[date_t]][q])

# Step 3: Create Basemap
m = b_map(0, -80, 360, 80, 'mill')
lons, lats = np.meshgrid(lons, lats)
m_map = Basemap(projection='mill', llcrnrlat=-90, urcrnrlat=90, llcrnrlon=0, urcrnrlon=360)
m_map.drawmapboundary(fill_color='LightCyan')
m_map.fillcontinents(color='LightSteelBlue',lake_color='LightCyan', zorder=0)

# Step 5: Plot Storm Track by day along with respective CO2 data from GOSAT satellite then save each day for each storm as a .png image file.
cb = 0
for st2 in range(len(storm_stats)-1):
    if len(storm_stats[str(st2)]['time']) > 3:
        if cb == 0: st2_00 = st2
        start_date = storm_stats[str(st2)]['time'][0]
        xlon, ylat = m(storm_stats[str(st2)]['lons'],storm_stats[str(st2)]['lats'])
        periods = len(storm_stats[str(st2)]['time'])
        days = len(storm_stats[str(st2)]['time']) / 4   # How many days did the storm last?
        if periods % 4 > 0 and periods > 4:
            days+=1 # Add an extra day if the storm isn't a multiple of 4 periods and > 1 Day
        zper=0  # Start at period 0
        startp = (((master_t[start_date]-data.variables['time'][0])/6)%4)
        for zday in range(days):    # Days loop
            txts = [False,False,False,False]
            if zper < periods:  # Only draw point if it exists
                if zper==0 and startp==0:
                    txt1 = plt.text(xlon[zper],ylat[zper],'*',fontsize=20,fontweight='heavy',color='k')
                    zper+=1         # Move on to the next period
                    txts[0]=True
                elif zper > 0:
                    txt1 = plt.text(xlon[zper],ylat[zper],'*',fontsize=20,fontweight='heavy',color='k')
                    zper+=1         # Move on to the next period
                    txts[0]=True
            if zper < periods:  # Only draw point if it exists
                if zper==0 and startp==1:
                    txt2 = plt.text(xlon[zper],ylat[zper],'*',fontsize=20,fontweight='heavy',color='k')
                    zper+=1         # Move on to the next period
                    txts[1]=True
                elif zper > 0:
                    txt2 = plt.text(xlon[zper],ylat[zper],'*',fontsize=20,fontweight='heavy',color='k')
                    zper+=1         # Move on to the next period
                    txts[1]=True
            if zper < periods:  # Only draw point if it exists
                if zper==0 and startp==2:
                    txt3 = plt.text(xlon[zper],ylat[zper],'*',fontsize=20,fontweight='heavy',color='k')
                    zper+=1         # Move on to the next period
                    txts[2]=True
                elif zper > 0:
                    txt3 = plt.text(xlon[zper],ylat[zper],'*',fontsize=20,fontweight='heavy',color='k')
                    zper+=1         # Move on to the next period
                    txts[2]=True
            if zper < periods:  # Only draw point if it exists
                if zper==0 and startp==3:
                    txt4 = plt.text(xlon[zper],ylat[zper],'*',fontsize=20,fontweight='heavy',color='k')
                    zper+=1         # Move on to the next period
                    txts[3]=True
                elif zper > 0:
                    txt4 = plt.text(xlon[zper],ylat[zper],'*',fontsize=20,fontweight='heavy',color='k')
                    zper+=1         # Move on to the next period
                    txts[3]=True
            
            date_string = str(hour_date(master_t,start_date+(zday*4)))
            co2_filename = '/home/meteo/cmm6461/meteo473/spring2015/milestone4/GOSAT_data/' + date_string[5:7] + '/' + date_string[8:10] + '.nc'
            
            # open OpenDAP dataset.
            co2data = Dataset(co2_filename)

            #co2 = co2data.variables['xCO2-bias_corrected'][:]
            co2 = co2data.variables['xCO2'][:]
            co2_lats = co2data.variables['latitude'][:]
            co2_lons = co2data.variables['longitude'][:]
            co2_time = co2data.variables['time'][:]
                
            master_co2 = []
            master_co2_lats = []
            master_co2_lons = []

            for f in range(len(co2)):
                co2[f]=co2[f]*(10**6)
                if co2_lons[f] < 0:
                    co2_lons[f] = co2_lons[f] +360
                if co2[f] > 380 and co2[f] < 410:
                    master_co2.append(co2[f])
                    master_co2_lats.append(co2_lats[f])
                    master_co2_lons.append(co2_lons[f])

            x,y = m_map(master_co2_lons, master_co2_lats)
            co2_data = plt.scatter(x,y, c=master_co2, s=20, cmap=plt.cm.jet,edgecolors='none')
            if cb == 0:
                plt.colorbar(pad=0.02, shrink=0.5)
                cb+=1
            plt.title(date_string[5:7] + '/' + date_string[8:10] + '_Storm' + str(st2) + '_Day' + str(zday))
            plt.savefig('strmco2_' + str(st2) + str(zday))
            if txts[0]: txt1.remove()
            if txts[1]: txt2.remove()
            if txts[2]: txt3.remove()
            if txts[3]: txt4.remove()
            co2_data.remove()

# STEP 6: Open webbrowser application to display image data.

slides = webbrowser.open('strmco2_' + str(st2_00) + '0.png','r+')


# STEP 7: Examine co2 data for 'x' number of days before/after a low pressure system is identified
#   Identify the mean of the co2 values for 'y' degrees around each storm point for each day
#   Identify the change in mean CO2 levels for 'x' number of days before/after a low pressure system is identified

# Loop to print out mean CO2 over storm track for 'co2_d' days before to 'co2_d' days after.
print "Storm " + strm_case + ':'
m_co2 = []
for ch_day in range(storm_stats[strm_case]['time'][0]-(co2_d*4),storm_stats[strm_case]['time'][0]+(co2_d*4)):
    mean_co2 = 0
    count = 0
    date_string = str(hour_date(master_t,ch_day))
    co2_filename = '/home/meteo/cmm6461/meteo473/spring2015/milestone4/GOSAT_data/' + date_string[5:7] + '/' + date_string[8:10] + '.nc'
    
    # open OpenDAP dataset.
    co2data = Dataset(co2_filename)

    #co2 = co2data.variables['xCO2-bias_corrected'][:]
    co2 = co2data.variables['xCO2'][:]
    co2_lats = co2data.variables['latitude'][:]
    co2_lons = co2data.variables['longitude'][:]
    co2_time = co2data.variables['time'][:]
        
    master_co2 = []
    master_co2_lats = []
    master_co2_lons = []

    for f in range(len(co2)):
        co2[f]=co2[f]*(10**6)
        if co2_lons[f] < 0:
            co2_lons[f] = co2_lons[f] +360
        if co2[f] > 380 and co2[f] < 410:
            master_co2.append(co2[f])
            master_co2_lats.append(co2_lats[f])
            master_co2_lons.append(co2_lons[f])

    storm = storm_stats[strm_case]
    for co2_lat in range(len(master_co2_lats)):
        for strm_lat in range(len(storm['lats'])):
            if abs(master_co2_lats[co2_lat]-(storm['lats'][strm_lat]) <= co2_area and abs(master_co2_lons[co2_lat])-storm['lons'][strm_lat]) <= co2_area:
                mean_co2 += master_co2[co2_lat]
                count += 1
    mean_co2 = mean_co2 / count
    m_co2.append(mean_co2)
    if len(m_co2) < 2 or m_co2[-1] != m_co2[-2]:    
        print date_string[5:7] + '/' + date_string[8:10] + ": " + str(mean_co2)
print "Change in mean co2 levels for " + str(co2_area) + " degrees of distance in all directions around the storm path is " + str(m_co2[-1]-m_co2[0]) + " from " + str(co2_d) + " days before/after the first/last recorded day of the system."
