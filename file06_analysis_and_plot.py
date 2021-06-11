#!/usr/bin/env python
# coding: utf-8

# In[46]:


get_ipython().run_line_magic('matplotlib', 'inline')
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib
import numpy as np
from numpy import logical_and as np_and
from numpy import logical_or as np_or
import HALO_processing as HALO
import copy
import csv
import math
import os
import time
import pandas as pd
import wind_calc
from scipy import stats
import pickle
from scipy.optimize import curve_fit
import uncertainties.unumpy as unp
import uncertainties as unc
import warnings
import matplotlib.font_manager as font_manager

font_dirs = ['fonts/', ]
font_files = font_manager.findSystemFonts(fontpaths=font_dirs)
font_list = font_manager.createFontList(font_files)
font_manager.fontManager.ttflist.extend(font_list)
plt.rcParams["font.family"] = "Arial"
# ignore warning
warnings.filterwarnings('ignore')


# In[47]:


def epa_data_input():
    epa_dict = pd.read_csv('Analysis_dataset/Taiwan_EPA_observation.csv')
    epa_dict = epa_dict[epa_dict.hour.isin([i for i in range(8, 18)])]
    epa_dict['CO'] = epa_dict['CO']*1000 #ppm to ppb
    epa_dict['SO2_over_O3'] = epa_dict['SO2']/epa_dict['O3']
    epa_dict['CO_over_O3'] = epa_dict['CO']/epa_dict['O3']
    epa_dict['NO2_over_O3'] = epa_dict['NO2']/epa_dict['O3']
    epa_dict['NO_over_O3'] = epa_dict['NO']/epa_dict['O3']
    epa_dict['CH4_over_O3'] = epa_dict['CH4']/epa_dict['O3']*1000 #ppm/ppb to ppb/ppb
    March_epa_dict = epa_dict[np_and(epa_dict.mm==3, epa_dict.dd.isin([i for i in range(1, 32)]))]
    epa_dict = epa_dict[np_or(np_and(epa_dict.mm==3, epa_dict.dd.isin([i for i in range(17, 32)])),
                                      np_and(epa_dict.mm==4, epa_dict.dd.isin([i for i in range(1, 8)])))]
    epa_loc_dict = pd.read_csv('Analysis_dataset/city_location.csv', encoding="Big5", index_col='Unnamed: 0')
    return epa_dict, March_epa_dict, epa_loc_dict


# In[48]:


def col_SO2O3_slope(epa_dict, col):
    place_list, slope_list, r2_list, p_value_list, slope_err_list = ([] for _ in range(5))
    slope_list_2, r2_list_2, p_value_list_2, slope_err_list_2 = ([] for _ in range(4))
    for p in set(epa_dict.station):
        temp_epa = epa_dict[epa_dict.station == p]
        try:
            mask_O3 = ~np.isnan(temp_epa['SO2_over_O3']) & ~np.isnan(temp_epa['%s_over_O3' %col])
            slope, intercept, r_value, p_value, std_err = stats.linregress(temp_epa['SO2_over_O3'][mask_O3], (temp_epa['%s_over_O3' %col])[mask_O3])
        except:
            slope, intercept, r_value, p_value, std_err = (np.nan for _ in range(5))
        try:
            mask = ~np.isnan(temp_epa['SO2']) & ~np.isnan(temp_epa['%s' %col])
            slope_2, intercept_2, r_value_2, p_value_2, std_err_2 = stats.linregress(temp_epa['SO2'][mask], (temp_epa['%s' %col])[mask])
        except:
            slope_2, intercept_2, r_value_2, p_value_2, std_err_2 = (np.nan for _ in range(5))
        place_list.append(p)
        slope_list.append(slope)
        r2_list.append(r_value**2)
        p_value_list.append(p_value)
        slope_err_list.append(std_err/np.sqrt(sum(mask_O3)))
        slope_list_2.append(slope_2)
        r2_list_2.append(r_value_2**2)
        p_value_list_2.append(p_value_2)
        slope_err_list_2.append(std_err_2/np.sqrt(sum(mask)))
    slope_output = pd.DataFrame({f'{col}/O3-SO2/O3_slope':slope_list,
                                 f'{col}/O3-SO2/O3_R2':r2_list,
                                 f'{col}/O3-SO2/O3_p_value':p_value_list,
                                 f'{col}/O3-SO2/O3_slope_err':slope_err_list,
                                 f'{col}-SO2_slope':slope_list_2,
                                 f'{col}-SO2_R2':r2_list_2,
                                 f'{col}-SO2_p_value':p_value_list_2,
                                 f'{col}-SO2_slope_err':slope_err_list_2,})
    return slope_output, place_list


# In[49]:


def slope_data_propcess(epa_dict):
    CH4_series, place_list = col_SO2O3_slope(epa_dict, 'CH4')
    CO_series, place_list = col_SO2O3_slope(epa_dict, 'CO')
    NO_series, place_list = col_SO2O3_slope(epa_dict, 'NO')
    NO2_series, place_list = col_SO2O3_slope(epa_dict, 'NO2')
    SO2_series, place_list = col_SO2O3_slope(epa_dict, 'SO2')
    slope_data = pd.concat([CH4_series, CO_series, NO_series, NO2_series, SO2_series], axis=1)
    slope_data.index = place_list
    slope_data['slope_x'] = 1/(slope_data['CO/O3-SO2/O3_slope']*0.16+                               slope_data['NO/O3-SO2/O3_slope']*24+                               slope_data['NO2/O3-SO2/O3_slope']*18.7+1)
    slope_data['slope_x_err'] = slope_data['slope_x']*np.sqrt((slope_data['CO/O3-SO2/O3_slope_err']*0.16)**2+                                                              (slope_data['NO/O3-SO2/O3_slope_err']*24)**2+                                                              (slope_data['NO2/O3-SO2/O3_slope_err']*18.7)**2)
    slope_data['p_value_all'] = np.logical_and(slope_data['CO/O3-SO2/O3_p_value']<0.05,
                                               np.logical_and(slope_data['NO/O3-SO2/O3_p_value']<0.05, slope_data['NO2/O3-SO2/O3_p_value']<0.05))
    return slope_data


# In[50]:


def pm_sinica_input():
    pm_sinica = pd.read_csv('Analysis_dataset/PM25_sinica_11cities.csv', encoding="Big5")
    pm_sinica = pm_sinica[pm_sinica['Time']=='D']
    return pm_sinica


# In[51]:


def aircraft_data_input():
    aircraft_data = pd.read_csv('Analysis_dataset/Western_Taiwan_dataset-1min.csv')
    aircraft_data['hour'] = aircraft_data['date'].astype(str) +' '+ aircraft_data['Local time'] 
    for header in aircraft_data.columns:
        if header not in ['date', 'Local time', 'place', 'hour']:
            aircraft_data[header] = pd.to_numeric(aircraft_data[header], errors='coerce')
        if header in ['SO2', 'NO', 'NOy', 'CO', 'CO2', 'O3_AMTEX', 'CH4', 'O3_FAIROCI']:
            aircraft_data.loc[aircraft_data[header] < 0] = np.nan
    aircraft_data['SO2'] /= 1000 # to ppb
    aircraft_data['NO2'] = aircraft_data['NOy']*(1-0.123) - aircraft_data['NO'] # minus PAN and NO
    aircraft_data['O3'] = np.nanmean([aircraft_data['O3_AMTEX'], aircraft_data['O3_FAIROCI']], axis=0)
    aircraft_data['SO2_over_CO'] = aircraft_data['SO2']/aircraft_data['CO']
    aircraft_data['SO2_over_O3'] = aircraft_data['SO2']/aircraft_data['O3']
    aircraft_data['CO_over_O3'] = aircraft_data['CO']/aircraft_data['O3']
    aircraft_data['hour'] = aircraft_data['Local time'].str[:2]
    aircraft_data['hour'] = aircraft_data['hour'].astype(float)
    aircraft_data = aircraft_data[np.logical_and(aircraft_data['hour']>=8, aircraft_data['hour']<=17)]
    aircraft_data = aircraft_data[aircraft_data['height_in_hundred']==600]
    aircraft_data['NO2_over_O3'] = aircraft_data['NO2']/aircraft_data['O3']
    aircraft_data['NO_over_O3'] = aircraft_data['NO']/aircraft_data['O3']
    aircraft_data['CH4_over_O3'] = aircraft_data['CH4']/aircraft_data['O3']*1000
    aircraft_data['date_str'] = aircraft_data['date'].apply(lambda x: str(int(x-20180000)))
    mask_NO_SO2_O3 = np.logical_and(aircraft_data['SO2_over_O3']>0.15, aircraft_data['NO2_over_O3']<0.1)
    aircraft_data['NO2_over_O3'][mask_NO_SO2_O3] = np.nan
    return aircraft_data


# In[52]:


def aircraft_all_input():
    date_set_all = ['20180308', '20180310', '20180312', '20180317', '20180319', '20180320',
                '20180322', '20180324', '20180326', '20180328', '20180330', '20180403',
                '20180404', '20180407', '20180409']

    date_hour_all = [(int(dd), h) for dd in date_set_all for h in range(8, 18)]

    aircraft_all = pd.read_csv('Analysis_dataset/All_area_dataset-15s.csv')
    aircraft_all['hour'] = aircraft_all['Local time'].str[:2]
    aircraft_all['hour'] = aircraft_all['hour'].astype(float)
    aircraft_all['date_str'] = aircraft_all['date'].apply(lambda x: str(int(x-20180000))).astype(str)
    for header in aircraft_all.columns:
        if header not in ['date', 'Local time', 'hour']:
            aircraft_all[header] = pd.to_numeric(aircraft_all[header], errors='coerce')
        if header in ['SO2', 'NO', 'NOy', 'CO', 'CO2', 'O3_AMTEX', 'CH4', 'O3_FAIROCI']:
            aircraft_all.loc[aircraft_all[header] < 0] = np.nan
    aircraft_all['SO2'] /= 1000 # to ppb
    #aircraft_all['NO2'] = aircraft_all['NOy']*(1-0.123) - aircraft_all['NO'] # minus PAN and NO
    aircraft_all['O3'] = np.nanmean([aircraft_all['O3_AMTEX'], aircraft_all['O3_FAIROCI']], axis=0)
    aircraft_all['SO2_over_O3'] = aircraft_all['SO2']/aircraft_all['O3']
    aircraft_all['CO_over_O3'] = aircraft_all['CO']/aircraft_all['O3']
    aircraft_all['NO_over_O3'] = aircraft_all['NO']/aircraft_all['O3']
    aircraft_all = aircraft_all[~np.isnan(aircraft_all.date)]
    return aircraft_all


# In[53]:


class Point():
    def __init__(self, lat, lon):
        self.lat = lat
        self.lon = lon
    def slope(self, p2):
        return (p2.lat-self.lat)/(p2.lon-self.lon)
    def intercept(self, p2):
        return (p2.lat-p2.lon*self.slope(p2))
class Line():
    def __init__(self, slope, intercept):
        self.slope = slope
        self.intercept = intercept
    def cal_lon(self, lat_position):
        return Point(lat_position, (lat_position-self.intercept)/self.slope)
            
def sub_region_test(lat, lon, region, boundary_only=False):
    diff = 0.35
    lat_list = np.array(sorted([diff*int(i) for i in range(6)], reverse=True)) + 22.8
    
    #test if the location is offshore
    p_ul = Point(24.72, 120.52)
    p_ll = Point(22.99, 119.85)
    p_ur = Point(24.60, 120.84)
    p_lr = Point(23.10, 120.28)
    
    L_ul  = Line(p_ul.slope(p_ll), p_ul.intercept(p_ll))
    L_ur = Line(p_ur.slope(p_lr), p_ur.intercept(p_lr))
    separation_L1 = Point(24.43, 120.56)
    separation_R1 = Point(24.42, 120.72)
    separation_L2 = Point(23.41, 120.02)
    separation_R2 = Point(23.33, 120.36)
    if region == 'region_a':
        p1 = L_ur.cal_lon(lat_list[1]) #pLR 
        p2 = L_ul.cal_lon(lat_list[1]) #pLL
        p3 = Point(lat_list[0], 120.54) #pUL
        p4 = Point(lat_list[0], 120.94) #pUR
    elif region == 'region_b':
        p1 = L_ur.cal_lon(lat_list[2]) #pLR 
        p2 = L_ul.cal_lon(lat_list[2]) #pLL
        p3 = L_ul.cal_lon(lat_list[1]) #pUL
        p4 = L_ur.cal_lon(lat_list[1]) #pUR
    elif region == 'region_c':
        p1 = L_ur.cal_lon(lat_list[3]) #pLR 
        p2 = L_ul.cal_lon(lat_list[3]) #pLL
        p3 = L_ul.cal_lon(lat_list[2]) #pUL
        p4 = L_ur.cal_lon(lat_list[2]) #pUR
    elif region == 'region_d':
        p1 = L_ur.cal_lon(lat_list[4]) #pLR 
        p2 = L_ul.cal_lon(lat_list[4]) #pLL
        p3 = L_ul.cal_lon(lat_list[3]) #pUL
        p4 = L_ur.cal_lon(lat_list[3]) #pUR
    elif region == 'region_e':
        p1 = Point(lat_list[5], 120.42) #pLR 
        p2 = Point(lat_list[5], 120.10) #pLL
        p3 = L_ul.cal_lon(lat_list[4]) #pUL
        p4 = L_ur.cal_lon(lat_list[4]) #pUR
        
    if boundary_only:
        return (p1, p2, p3, p4)
    else:
        L1 = Line(p1.slope(p2), p1.intercept(p2))
        L2 = Line(p2.slope(p3), p2.intercept(p3))
        L3 = Line(p3.slope(p4), p3.intercept(p4))
        L4 = Line(p4.slope(p1), p4.intercept(p1))
        if L1.slope >= 0:
            cal_L1 = (L1.slope*lon - lat + L1.intercept) <= 0
        else:
            cal_L1 = (L1.slope*lon - lat + L1.intercept) <= 0
        if L2.slope > 0:
            cal_L2 = (L2.slope*lon - lat + L2.intercept) >= 0
        else:
            cal_L2 = (L2.slope*lon - lat + L2.intercept) <= 0
        if L3.slope <= 0:
            cal_L3 = (L3.slope*lon - lat + L3.intercept) >= 0
        else:
            cal_L3 = (L3.slope*lon - lat + L3.intercept) >= 0
        if L4.slope > 0:
            cal_L4 = (L4.slope*lon - lat + L4.intercept) <= 0
        else:
            cal_L4 = (L4.slope*lon - lat + L4.intercept) >= 0
        return np.logical_and(np.logical_and(cal_L1, cal_L2), np.logical_and(cal_L3, cal_L4))


# In[54]:


def func(x, a, b):
    """The fitting function"""
    return a*(np.array(x))+b


# In[55]:


def map_region_border(p1, p2, p3, p4, ax):
    line_setting_2 = dict(linestyle='--', linewidth=2.5)
    ax.plot([p1.lon, p2.lon],[p1.lat, p2.lat], **line_setting_2, color='r')
    ax.plot([p2.lon, p3.lon],[p2.lat, p3.lat], **line_setting_2, color='r')
    ax.plot([p3.lon, p4.lon],[p3.lat, p4.lat], **line_setting_2, color='r')
    ax.plot([p4.lon, p1.lon],[p4.lat, p1.lat], **line_setting_2, color='r')
    return None


# In[56]:


def epa_col_to_SO2O3_slope_cal(epa_input, col, slope_list, r2_list):
    try:
        mask = ~np.isnan(epa_input['SO2_over_O3']) & ~np.isnan(epa_input[col])
        slope, intercept, r_value, p_value, std_err = stats.linregress(epa_input['SO2_over_O3'][mask],
                                                                       epa_input[col][mask])
        if slope > 0 and p_value < 0.05:
            slope_list.append(slope)
        else:
            slope_list.append(np.nan)
        r2_list.append(r_value**2)
    except:
        slope_list.append(np.nan)
        r2_list.append(np.nan)


# In[57]:


def air_col_to_SO2O3_slope_cal(air_input, col, slope_list, slope_err_list):
    mask = ~np.isnan(air_input['SO2_over_O3']) & ~np.isnan(air_input[col])
    slope, intercept, r_value, p_value, std_err = stats.linregress(air_input['SO2_over_O3'][mask],
                                                                   air_input[col][mask])
    if p_value < 0.05:
        slope_list.append(slope)
        slope_err_list.append(std_err)
    else:
        slope_list.append(np.nan)
        slope_err_list.append(np.nan)
    return sum(mask)


# In[58]:


def plot_figure_2(epa_dict):
    """
    Plot Figure 2: The CO-SO2 (a, b), NO-SO2 (c, d), and NO2-SO2 (e, f) relationships 
    before (a, c, e) and after (b, d, f) the O3 normalization of Tainan station.
    """
    fig, ((ax00, ax01), (ax10, ax11), (ax20, ax21) ) =             plt.subplots(3, 2, sharex=False, sharey=False, figsize=(14, 16))
    fig.tight_layout(pad=6.0)

    mark_size = 12
    xytick_size = 16
    label_size = 22
    scatter_arg = dict(s=25, c='tab:blue', alpha=0.5, marker='o')
    p = '臺南'

    index = 0
    ax_list = [ax00, ax01, ax10, ax11, ax20, ax21]
    plot_vars = [('SO2', 'CO'), ('SO2_over_O3', 'CO_over_O3'),
                 ('SO2', 'NO'), ('SO2_over_O3', 'NO_over_O3'),
                 ('SO2', 'NO2'), ('SO2_over_O3', 'NO2_over_O3')]
    x_max = [12, 0.3]*3
    y_max = [1000, 40, 15, 0.8, 35, 1.5]
    y_label_list = ['CO (ppb)', 'CO/$\mathregular{O_3}$ ratio (ppb/ppb)',
                    'NO (ppb)', 'NO/$\mathregular{O_3}$ ratio (ppb/ppb)',
                    '$\mathregular{NO_2}$ (ppb)', '$\mathregular{NO_2}$/$\mathregular{O_3}$ ratio (ppb/ppb)']
    label_text = ['a', 'b', 'c', 'd', 'e', 'f']
    for index in range(6):
        (x_col, y_col) = plot_vars[index]
        ax = ax_list[index]
        ax.set_ylim(0, y_max[index])
        ax.set_xlim(0, x_max[index])
        (xmin, xmax), (ymin, ymax) = ax.get_xlim(), ax.get_ylim()
        if x_col == 'SO2':
            ax.set_xlabel('$\mathregular{SO_2}$ (ppb)', color='k', fontsize=label_size, labelpad=10)
        else:
            ax.set_xlabel('$\mathregular{SO_2}$/$\mathregular{O_3}$ ratio (ppb/ppb)', color='k', fontsize=label_size, labelpad=10)
        ax.set_ylabel(y_label_list[index], color='k', fontsize=label_size, labelpad=10)
        temp_epa = epa_dict[epa_dict.station == p]
        ax.scatter(temp_epa[x_col], temp_epa[y_col], **scatter_arg)
        mask = ~np.isnan(temp_epa[x_col]) & ~np.isnan(temp_epa[y_col])
        slope_2, intercept_2, r_value_2, p_value_2, std_err_2 = stats.linregress(temp_epa[x_col][mask], temp_epa[y_col][mask])
        y_low, y_high = xmin*slope_2+intercept_2, xmax*slope_2+intercept_2
        ax.plot([xmin, xmax], [y_low, y_high], c='r', linewidth=2, zorder=1)           
        ax.text((xmax-xmin)*0.03+xmin, (ymax-ymin)*0.9+ymin, f'({label_text[index]})', fontsize=label_size, weight='bold')
        ax.text((xmax-xmin)*0.75+xmin, (ymax-ymin)*0.9+ymin, '$\mathregular{R^2}$: %.3f' %r_value_2**2, fontsize=label_size-3)
        ax.xaxis.set_tick_params(labelsize=xytick_size)
        ax.yaxis.set_tick_params(labelsize=xytick_size)    
        index += 1
    fig.savefig('figures/Fig.2-before_after_division_O3.jpg', dpi=300, bbox_inches='tight')
    plt.show()


# In[59]:


def plot_sfc_date_range_avg_sulfate(ddd, ax, label_text, 
                                    pm_sinica, March_epa_dict, p_to_en,
                                    legend=False, detail=False):
    label_fontsize = 26
    legend_dict = dict(ncol=2, bbox_to_anchor=(0.5, -0.22), loc='center', borderaxespad=0, labelspacing=0.5, frameon=True, fontsize=28)
    epa_scatter_s = 100
    epa_scatter_alpha = 0.5
    ytick_size=18
    textsize = 18

    place_list, pm_list, pm_err, r2_list = ([] for _ in range(4))
    slope_COO3_list, slope_NOO3_list, slope_NO2O3_list = ([] for _ in range(3))
    r2_COO3_list, r2_NOO3_list, r2_NO2O3_list = ([] for _ in range(3))
    SO2_list, RH_list = ([] for _ in range(2))

    for p in set(pm_sinica.place):
        temp_epa = March_epa_dict[np_and(March_epa_dict.station == p, March_epa_dict.date.isin(ddd))]
        place_list.append(p)
        SO2_list.extend(temp_epa.SO2)
        RH_list.extend(temp_epa.RH)
        pm_sinica_mask = (pm_sinica['Date'].isin(ddd)) & (pm_sinica['place']==p)
        pm_list.append(float(pm_sinica[pm_sinica_mask]['SO4'].mean()))
        pm_err.append(pm_sinica[pm_sinica_mask]['SO4'].std()/np.sqrt(pm_sinica[pm_sinica_mask]['SO4'].count()))
        epa_col_to_SO2O3_slope_cal(temp_epa, 'CO_over_O3', slope_COO3_list, r2_COO3_list)
        epa_col_to_SO2O3_slope_cal(temp_epa, 'NO_over_O3', slope_NOO3_list, r2_NOO3_list)
        epa_col_to_SO2O3_slope_cal(temp_epa, 'NO2_over_O3', slope_NO2O3_list, r2_NO2O3_list)
    pm_data = pd.DataFrame({'slope_COO3':slope_COO3_list, 'r2_COO3':r2_COO3_list,
                            'slope_NOO3':slope_NOO3_list, 'r2_NOO3':r2_NOO3_list,
                            'slope_NO2O3':slope_NO2O3_list, 'r2_NO2O3':r2_NO2O3_list,
                            'pm':pm_list, 'pm_err':pm_err, 'place':place_list})

    pm_data['slope_x'] = 1/(pm_data['slope_COO3']*0.16+pm_data['slope_NOO3']*24+pm_data['slope_NO2O3']*18.7+1)
    c_l = [plt.cm.rainbow(i/11) for i in range(11)]
    place_list = list(sorted(set(place_list)))
    place_order = [pm_data.index[pm_data['place'] == p][0] for p in place_list]
    
    for i in range(len(set(place_list))):
        select = pm_data['place'] == place_list[i]
        ax.scatter(pm_data.loc[select, 'slope_x'], pm_data.loc[select, 'pm'],
                   s=50, color=c_l[i], alpha=1, marker='o', 
                   label=p_to_en[place_list[i]], zorder=5)
    if legend:
        ax.legend(loc='center left', bbox_to_anchor=(1.05, 0.5), fontsize=16)
        
    ax.errorbar(pm_data.loc[place_order, 'slope_x'], pm_data.loc[place_order, 'pm'],
                yerr=pm_data.loc[place_order, 'pm_err'],
                linestyle="None", fmt='-o', ecolor='darkgrey', c='k', markersize=8, capsize=5, zorder=3)
    
    mask = ~np.isnan(pm_data['slope_x']) & ~np.isnan(pm_data['pm'])
    slope_final, intercept_final, r_value_final, p_value_final, std_err_final = stats.linregress(pm_data['slope_x'][mask], pm_data['pm'][mask])
    ax.set_ylabel('Sulfate ($\mathregular{{\mu}g}$/$\mathregular{m^3}$)', color='k', fontsize = label_fontsize, labelpad=10)
    ax.set_xlabel('$\mathregular{SO_2}$ oxidation rate fraction', color='k', fontsize = label_fontsize, labelpad=10)
    ax.xaxis.set_tick_params(labelsize=ytick_size)
    ax.yaxis.set_tick_params(labelsize=ytick_size)
        
    n = len(pm_data['pm'][mask])
    popt, pcov = curve_fit(func, pm_data['slope_x'][mask], pm_data['pm'][mask])

    # retrieve parameter values to compute R2
    a, b = popt[0], popt[1]
    r2 = 1.0-(sum((pm_data['pm'][mask]-func(pm_data['slope_x'][mask],a,b))**2)/((n-1.0)*np.var(pm_data['pm'][mask], ddof=1)))
    # calculate parameter confidence interval
    a, b = unc.correlated_values(popt, pcov)

    
    # calculate explicable SO4 percentage
    SO4_max = pm_data['pm'].max()
    SO4_max_err = pm_data[pm_data['pm'] == pm_data['pm'].max()]['pm_err']
    explicable_SO4_percentage = (SO4_max-intercept_final)/SO4_max*100
    rxn_time = (SO4_max-intercept_final)/(1.5e-12*8.3e6*2.46e10*3600/6.02e23*1e6*96*1e6)/np.nanmean(SO2_list)
    
    # calculate regression confidence interval
    xmin, xmax, (ymin, ymax) = pm_data['slope_x'].min(), pm_data['slope_x'].max(), ax.get_ylim()
    ax.set_xlim(xmin/1.15, xmax*1.15)
    ax.set_ylim(ymin, ymax+0.5)
    px = np.linspace(xmin-0.005, xmax+0.005, num=50, endpoint=True)
    py = a*px+b
    nom = unp.nominal_values(py)
    std = unp.std_devs(py)
    
    # plot the regression line and uncertainty band (95% confidence)
    ax.plot(px, nom, c='r')
    ax.fill_between(px, nom - 1.96 * std, nom + 1.96 * std, color='orange', alpha=0.2)

    (_, xmax), (ymin, ymax) = ax.get_xlim(), ax.get_ylim()
    xmin = 0
    ax.set_xlim(0, xmax)
    ax.text((xmax-xmin)*0.55+xmin, (ymax-ymin)*0.23+ymin,  r'$\mathregular{slope:}$', fontsize=textsize)
    ax.text((xmax-xmin)*0.77+xmin, (ymax-ymin)*0.23+ymin,  r'$\mathregular{%.1f \pm %.1f}$' %(slope_final, pcov[0, 0]**0.5), fontsize=textsize)
    ax.text((xmax-xmin)*0.55+xmin, (ymax-ymin)*0.17+ymin, r'$\mathregular{intercept:}$', fontsize=textsize)
    ax.text((xmax-xmin)*0.77+xmin, (ymax-ymin)*0.17+ymin, r'$\mathregular{%.2f \pm %.2f}$' %(intercept_final, pcov[1, 1]**0.5), fontsize=textsize)
    ax.text((xmax-xmin)*0.55+xmin, (ymax-ymin)*0.09+ymin, r'$\mathregular{R^2:}$', fontsize=textsize)
    ax.text((xmax-xmin)*0.77+xmin, (ymax-ymin)*0.09+ymin, r'$\mathregular{%.3f}$' %(r_value_final**2), fontsize=textsize)
    ax.text((xmax-xmin)*0.55+xmin, (ymax-ymin)*0.03+ymin, '$\mathregular{p-value:}$', fontsize=textsize)  
    if p_value_final > 0.001:
        ax.text((xmax-xmin)*0.77+xmin, (ymax-ymin)*0.03+ymin, r'$\mathregular{%.3f}$' %(p_value_final), fontsize=textsize)
    else:
        ax.text((xmax-xmin)*0.77+xmin, (ymax-ymin)*0.03+ymin, r'$\mathregular{%.1e}$' %(p_value_final), fontsize=textsize)
    ax.text((xmax-xmin)*0.03+xmin, (ymax-ymin)*0.91+ymin, f'{label_text}',
            fontsize=label_fontsize, weight='bold') 
    
    if detail:
        print('-'*10+'Surface'+'-'*10)
        print('Period:', ddd)
        print('slope: %.1f+/-%.1f,' %(slope_final, pcov[0, 0]**0.5))
        print('intercept: %.2f+/-%.2f,' %(intercept_final, pcov[1, 1]**0.5))
        print('R2:%.3f, p value: %.3E' %(r_value_final**2, p_value_final))
        print('explicable SO4 (ug/m3): %.1f' %(SO4_max-intercept_final))
        print('explicable SO4 percentage: %.1f+/-%.1f' %(explicable_SO4_percentage, pcov[1, 1]**0.5/SO4_max*100))
        print('RH: %.1f +/- %.1f' %(np.nanmean(RH_list), np.nanstd(RH_list)/np.sqrt(len(RH_list))))
        print('SO2: %.2f +/- %.2f' %(np.nanmean(SO2_list), np.nanstd(SO2_list)/np.sqrt(len(SO2_list))))
        print('Estimated reaction time: %.1f h' %rxn_time)
        print('Max rate fraction:', pm_data['slope_x'].max())


# In[60]:


def plot_airborne_avg_sulfate_to_so2_oxi_rate_frac(ax, label_text, aircraft_data, detail=False):
    label_fontsize = 26
    slope_CO_list, slope_CO_err_list,     slope_NO_list, slope_NO_err_list,     slope_NO2_list, slope_NO2_err_list = ([] for _ in range(6))
    so4, so4_err, SO2_list, RH_list = ([] for _ in range(4))

    for region_name in ['region_a', 'region_b', 'region_c', 'region_d', 'region_e']:
        temp_aircraft = aircraft_data[sub_region_test(aircraft_data.lat, aircraft_data.lon, region=region_name)]
        so4.append(temp_aircraft['SO4'].mean())
        so4_err.append(temp_aircraft['SO4'].std()/np.sqrt(temp_aircraft['SO4'].count()))
        SO2_list.extend(temp_aircraft['SO2'])
        RH_list.extend(temp_aircraft['RELHUM'])
        count_CO = air_col_to_SO2O3_slope_cal(temp_aircraft, 'CO_over_O3', slope_CO_list, slope_CO_err_list)
        count_NO = air_col_to_SO2O3_slope_cal(temp_aircraft, 'NO_over_O3', slope_NO_list, slope_NO_err_list)
        count_NO2 = air_col_to_SO2O3_slope_cal(temp_aircraft, 'NO2_over_O3', slope_NO2_list, slope_NO2_err_list)    
    slope_x = 1/(np.array(slope_CO_list)*0.16+np.array(slope_NO_list)*24+np.array(slope_NO2_list)*18.7+1)  
    slope_x_err = slope_x*np.sqrt((0.16*np.array(slope_CO_err_list)/np.sqrt(count_CO))**2+                                  (24*np.array(slope_NO_err_list)/np.sqrt(count_NO))**2+                                  (18.7*np.array(slope_NO2_err_list)/np.sqrt(count_NO2))**2)
    # -------------- Start Plotting --------------------
    ax.errorbar(slope_x, so4, yerr=so4_err, linestyle="None", fmt='-o', ecolor='darkgrey', c='k', markersize=8, capsize=5, zorder=3)
    c_l = [plt.cm.rainbow(i/4) for i in range(5)]
    region_list = ['region_a', 'region_b', 'region_c', 'region_d', 'region_e']
    for i in range(5):
        ax.scatter(slope_x[i], so4[i], s=50, color=c_l[i], alpha=1, marker='o', 
                   label=' '.join(['Region', region_list[i][-1].capitalize()]), zorder=5)
    ax.legend(loc='center left', bbox_to_anchor=(1.05, 0.5), fontsize=16)
    ax.xaxis.set_tick_params(labelsize=18)
    ax.yaxis.set_tick_params(labelsize=18)
    ax.set_xlabel('$\mathregular{SO_2}$ oxidation rate fraction', color='k', fontsize=22, labelpad=10)
    ax.set_ylabel('Sulfate ($\mathregular{{\mu}g}$/$\mathregular{m^3}$)', color='k', fontsize=22, labelpad=10)
    
    slope_final, intercept_final, r_value_final, p_value_final, std_err_final = stats.linregress(slope_x, so4)
    x = np.linspace(slope_x.min(), slope_x.max(), num=50, endpoint=True)
    n = len(so4)
    popt, pcov = curve_fit(func, slope_x, so4)

    # retrieve parameter values to compute R2
    a, b = popt[0], popt[1]
    r2 = 1.0-(sum((so4-func(slope_x,a,b))**2)/((n-1.0)*np.var(so4, ddof=1)))
    
    # calculate parameter confidence interval
    a, b = unc.correlated_values(popt, pcov)

    # calculate explicable SO4 percentage
    SO4_max = np.array(so4).max()
    SO4_max_err = np.array(so4).std()
    explicable_SO4_percentage = (SO4_max-intercept_final)/SO4_max*100
    rxn_time = (SO4_max-intercept_final)/(1.5e-12*8.3e6*2.46e10*3600/6.02e23*1e6*96*1e6)/np.nanmean(SO2_list)

    # calculate regression confidence interval
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    ax.set_xlim(xmin-0.005, xmax+0.005)
    ax.set_ylim(ymin-0.1, ymax+0.2)
    px = np.linspace(xmin-0.005, xmax+0.005, num=50, endpoint=True)
    py = a*px+b
    nom = unp.nominal_values(py)
    std = unp.std_devs(py)

    # plot the regression line and uncertainty band (95% confidence)
    ax.plot(px, nom, c='r')
    ax.fill_between(px, nom - 1.96 * std, nom + 1.96 * std, color='orange', alpha=0.2)

    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    textsize = 18
    ax.text((xmax-xmin)*0.56+xmin, (ymax-ymin)*0.23+ymin,  r'$\mathregular{slope:}$', fontsize=textsize)
    ax.text((xmax-xmin)*0.78+xmin, (ymax-ymin)*0.23+ymin,  r'$\mathregular{%.1f \pm %.1f}$' %(slope_final, pcov[0, 0]**0.5), fontsize=textsize)
    ax.text((xmax-xmin)*0.56+xmin, (ymax-ymin)*0.17+ymin, r'$\mathregular{intercept:}$', fontsize=textsize)
    ax.text((xmax-xmin)*0.78+xmin, (ymax-ymin)*0.17+ymin, r'$\mathregular{%.2f \pm %.2f}$' %(intercept_final, pcov[1, 1]**0.5), fontsize=textsize)
    ax.text((xmax-xmin)*0.56+xmin, (ymax-ymin)*0.09+ymin, r'$\mathregular{R^2:}$', fontsize=textsize)
    ax.text((xmax-xmin)*0.78+xmin, (ymax-ymin)*0.09+ymin, r'$\mathregular{%.3f}$' %(r_value_final**2), fontsize=textsize)
    ax.text((xmax-xmin)*0.56+xmin, (ymax-ymin)*0.03+ymin, '$\mathregular{p-value:}$', fontsize=textsize)  
    if p_value_final > 0.001:
        ax.text((xmax-xmin)*0.78+xmin, (ymax-ymin)*0.03+ymin, r'$\mathregular{%.3f}$' %(p_value_final), fontsize=textsize)
    else:
        ax.text((xmax-xmin)*0.78+xmin, (ymax-ymin)*0.03+ymin, r'$\mathregular{%.1e}$' %(p_value_final), fontsize=textsize)
    ax.text((xmax-xmin)*0.03+xmin, (ymax-ymin)*0.91+ymin, f'{label_text}',
            fontsize=label_fontsize, weight='bold') 
    
    if detail:
        print('-'*10+'Aircraft'+'-'*10)
        print('slope: %.1f+/-%.1f,' %(slope_final, pcov[0, 0]**0.5))
        print('intercept: %.2f+/-%.2f,' %(intercept_final, pcov[1, 1]**0.5))
        print('R2:%.3f, p value: %.3E' %(r_value_final**2, p_value_final))
        print('explicable SO4: %.2f' %(SO4_max-intercept_final))                             
        print('explicable SO4 percentage: %.1f+/-%.1f' %(explicable_SO4_percentage, pcov[1, 1]**0.5/SO4_max*100))
        print('RH: %.1f +/- %.1f' %(np.nanmean(RH_list), np.nanstd(RH_list)/np.sqrt(len(RH_list))))
        print('SO2: %.2f +/- %.2f' %(np.nanmean(SO2_list), np.nanstd(SO2_list)/np.sqrt(len(SO2_list))))
        print('Estimated reaction time: %.1f h' %rxn_time)
        print('Max rate fraction:', slope_x.max())


# In[61]:


def plot_figure_3(pm_sinica, March_epa_dict, p_to_en, aircraft_data, detail=False):
    """
    Figure 3: Relationships between sulfate concentrations and SO2 oxidation rate fraction
    of surface EPA stations during March 13th to 31st, 2018 (a) and airborne observation
    from March 17th to April 7th, 2018 (b). Fitting lines are plotted in solid red lines
    with orange shades of 95% confidential interval.
    """
    fig = plt.figure(figsize=(8, 12))
    ax1 = fig.add_axes([0.1, 0.1/2+1/2, 0.8, 0.8/2])
    ax2 = fig.add_axes([0.1, 0.1/2, 0.8, 0.8/2])

    # Surface period: '3/13-31'
    date_all = set(pm_sinica['Date'])
    plot_sfc_date_range_avg_sulfate(date_all, ax1, '(a)', pm_sinica, March_epa_dict, p_to_en, legend=True, detail=detail)
    # Airborne
    plot_airborne_avg_sulfate_to_so2_oxi_rate_frac(ax2, '(b)', aircraft_data, detail=detail)

    # save figure
    fig.savefig('figures/Fig.3-sulfate_SO4_oxi_rate_frac_sfc_and_air_western_Taiwan.jpg', dpi=300, bbox_inches='tight')
    plt.show()


# In[62]:


def plot_figure_4(pm_sinica, March_epa_dict, p_to_en):
    """
    Figure 4: The relationships between sulfate concentration and SO2 oxidation rate fraction
    using equation 9 based on surface measurements. The durations of (a), (b), and (c) are
    March 13rd – 19th, Mach 21st - 23rd, and Mach 25th – 28th, respectively.
    """
    fig = plt.figure(figsize=(8, 18))
    ax3 = fig.add_axes([0.1, 0.1/3, 0.8, 0.8/3])
    ax2 = fig.add_axes([0.1, 0.1/3+1/3, 0.8, 0.8/3])
    ax1 = fig.add_axes([0.1, 0.1/3+2/3, 0.8, 0.8/3])

    # period: '3/13-19'
    date_a = ['20180313', '20180314', '20180315', '20180316', '20180317', '20180318', '20180319']
    plot_sfc_date_range_avg_sulfate(date_a, ax1, '(a)', pm_sinica, March_epa_dict, p_to_en)
    # period: '3/21-23'
    date_b = ['20180321', '20180322', '20180323']
    plot_sfc_date_range_avg_sulfate(date_b, ax2, '(b)', pm_sinica, March_epa_dict, p_to_en)
    # period: '3/25-28'
    date_c = ['20180325', '20180326', '20180327', '20180328']
    plot_sfc_date_range_avg_sulfate(date_c, ax3, '(c)', pm_sinica, March_epa_dict, p_to_en)
    
    # save figure
    fig.savefig('figures/Fig.4-sulfate_SO4_oxi_rate_fraction_sfc_period.jpg', dpi=300, bbox_inches='tight')
    plt.show()


# In[63]:


def map_with_parallel_meridian(ax1, lon_min, lon_max, lat_min, lat_max, lonlat_gap):
    m = Basemap(projection='cyl',resolution='h',llcrnrlat=lat_min,         urcrnrlat=lat_max, llcrnrlon=lon_min, urcrnrlon=lon_max, ax=ax1)
    m.arcgisimage(service='World_Shaded_Relief', xpixels=1000, verbose=False, alpha=0.1)
    # draw parallels and meridians.
    parallels = np.arange(lat_min, lat_max+0.1, 2)
    m.drawparallels(parallels, labels=[False,False,False,False], fontsize=20)
    meridians = np.arange(lon_min, lon_max+0.1, 2)
    m.drawmeridians(meridians, labels=[False,False,False,False], fontsize=20)
    parallels_2 = np.arange(lat_min, lat_max+0.1, lonlat_gap)
    m.drawparallels(parallels_2, labels=[False,False,False,False], linewidth=0.5)
    meridians_2 = np.arange(lon_min, lon_max+0.1, lonlat_gap)
    m.drawmeridians(meridians_2, labels=[False,False,False,False], linewidth=0.5)
    return m

def ax_setting(ax, xmin, xmax, xinterval):
    ax.xaxis.set_tick_params(labelsize=16)
    ax.yaxis.set_tick_params(labelsize=16)
    ax.set_xlim(xmin, xmax+xinterval//2)


# In[64]:


def height_slope_package(aircraft_all, ax_left, ax_right, date, height, xmin, xmax, xinterval,
                         box_w, lon_min, lon_max, lon_tick_num, lat_min, lat_max, lat_tick_num,
                         ymin, ymax, lonlat_gap, label_left, label_right, detail=False):
    temp_dd = aircraft_all[np_and(aircraft_all['date']==date,  aircraft_all['height_in_hundred']==height)]
    temp_dd = temp_dd[np_and(temp_dd['lon']>=lon_min, temp_dd['lon']<=lon_max)]
    temp_dd = temp_dd[np_and(temp_dd['lat']>=lat_min, temp_dd['lat']<=lat_max)]
    
    # plot flight track (left panel)
    cmap = matplotlib.cm.get_cmap('rainbow')
    m = map_with_parallel_meridian(ax_left, lon_min, lon_max, lat_min, lat_max, lonlat_gap)
    x_interval = np.linspace(lon_min, lon_max, num=lon_tick_num,)
    y_interval = np.linspace(lat_min, lat_max, num=lat_tick_num,)
    ax_left.set_xticks(x_interval)
    ax_left.set_yticks(y_interval)
    ax_left.set_xticklabels(['%.1f$\degree$E' %i for i in x_interval])
    ax_left.set_yticklabels(['%.1f$\degree$N' %i for i in y_interval])
    ax_left.xaxis.set_tick_params(labelsize=20, pad=10)
    ax_left.yaxis.set_tick_params(labelsize=20, pad=10)
    index = ax_left.text(0.05, 0.95, label_left, transform=ax_left.transAxes,
                         fontsize=22, fontweight='bold', va='top')
    index.set_bbox(dict(facecolor='white', alpha=0.5, edgecolor='white'))  
    flight = m.scatter(temp_dd['lon'], temp_dd['lat'], c='r', label='Flights', s=5)

    # plot sulfate-SO2 oxidation rate fraction (right panel)
    # NO2O3-SO2O3
    mask = ~np.isnan(temp_dd['NOy']) & ~np.isnan(temp_dd['PAN_PPB'])
    slope_PAN, intercept_PAN, r_value_PAN, p_value_PAN, std_err_PAN = stats.linregress(temp_dd['NOy'][mask],
                                                                                       temp_dd['PAN_PPB'][mask])
    #temp_aircraft['NO2'] = temp_aircraft['NOy']*(1-0.123) - temp_aircraft['NO']
    temp_dd['NO2'] = temp_dd['NOy']*(1-slope_PAN) - temp_dd['NO']
    temp_dd['NO2_over_O3'] = temp_dd['NO2']/temp_dd['O3']
    mask_NO_SO2_O3 = np_and(temp_dd['SO2_over_O3']>0.15, temp_dd['NO2_over_O3']<0.1)
    temp_dd['NO2_over_O3'][mask_NO_SO2_O3] = np.nan
    
    lon_left = np.arange(lon_min, lon_max, lonlat_gap)
    lat_bottom = np.arange(lat_min, lat_max, lonlat_gap)
    so4_series, so4_mean, so4_err, so4_median = [], [], [], []
    slope_CO_list, slope_err_CO_list, slope_NO_list, slope_err_NO_list, slope_NO2_list, slope_err_NO2_list = ([] for _ in range(6))
    so2_all, RH = [], []
    
    for i in lon_left:
        for j in lat_bottom:
            temp = temp_dd[(temp_dd['lat'] > j) & (temp_dd['lat'] <= j+lonlat_gap)]
            if len(temp) > 0:
                temp = temp[(temp['lon'] > i) & (temp['lon'] <= i+lonlat_gap)]
                mask_2 = ~np.isnan(temp['SO2_over_O3']) & ~np.isnan(temp['CO_over_O3']) &                         ~np.isnan(temp['NO_over_O3']) & ~np.isnan(temp['NO2_over_O3'])
                num = mask_2.sum()
                if num > 2:
                    slope_CO, intercept_CO, r_value_CO, p_value_CO, std_err_CO = stats.linregress(temp['SO2_over_O3'][mask_2], temp['CO_over_O3'][mask_2])
                    slope_NO, intercept_NO, r_value_NO, p_value_NO, std_err_NO = stats.linregress(temp['SO2_over_O3'][mask_2], temp['NO_over_O3'][mask_2])
                    slope_NO2, intercept_NO2, r_value_NO2, p_value_NO2, std_err_NO2 = stats.linregress(temp['SO2_over_O3'][mask_2], temp['NO2_over_O3'][mask_2])
                    if  slope_CO>0 and slope_NO>0 and slope_NO2>0:
                        slope_CO_list.append(slope_CO)
                        slope_err_CO_list.append(std_err_CO)
                        slope_NO_list.append(slope_NO)
                        slope_err_NO_list.append(std_err_NO)
                        slope_NO2_list.append(slope_NO2)
                        slope_err_NO2_list.append(std_err_NO2)
                        so4_series.append(list(temp['SO4'].dropna()))
                        so4_mean.append(temp['SO4'].mean())
                        so4_median.append(temp['SO4'].median())
                        so4_series.append(list(temp['SO4'].dropna()))
                        so4_err.append(temp['SO4'].std()/np.sqrt(num))
                        so2_all.extend(temp['SO2'].dropna())
                        RH.extend(temp['RELHUM'].dropna())
    so4 = np.array(so4_mean)
    so4_err = np.array(so4_err)
    slope_x = 1/(np.array(slope_CO_list)*0.16+np.array(slope_NO_list)*24+np.array(slope_NO2_list)*18.7+1)
    mask_3 = ~np.isnan(slope_x) & ~np.isnan(so4)
    slope_3, intercept_3, r_value_3, p_value_3, std_err_3 = stats.linregress(slope_x[mask_3], so4[mask_3])

    ax_setting(ax_right, xmin, xmax, xinterval)
    n = len(so4[mask_3])
    popt, pcov = curve_fit(func, slope_x[mask_3], so4[mask_3])
    # retrieve parameter values to compute R2
    a, b = popt[0], popt[1]
    r2 = 1.0-(sum((so4[mask_3]-func(slope_x[mask_3],a,b))**2)/((n-1.0)*np.var(so4[mask_3],ddof=1)))
    # calculate parameter confidence interval
    a, b = unc.correlated_values(popt, pcov)

    # calculate explicable SO4 percentage
    SO4_max = so4[mask_3].max()
    SO4_max_err = so4[mask_3].std()
    explicable_SO4_percentage = (SO4_max-intercept_3)/SO4_max*100     
    rxn_time = (SO4_max-intercept_3)/(1.5e-12*8.3e6*2.46e10*3600/6.02e23*1e6*96*1e6)/np.nanmean(so2_all)
    
    # plot data
    ax_right.errorbar(slope_x, so4, yerr=so4_err, linestyle="None", fmt='-o', ecolor='darkgrey', c='k', markersize=8, capsize=5)
    # calculate regression confidence interval
    px = np.linspace(0, slope_x[mask_3].max()+0.1, num=50, endpoint=True) 
    py = a*px+b
    nom = unp.nominal_values(py)
    std = unp.std_devs(py)
    # plot the regression line and uncertainty band (95% confidence)
    ax_right.plot(px, nom, c='r')
    ax_right.fill_between(px, nom - 1.96 * std, nom + 1.96 * std, color='orange', alpha=0.2)
    # text setting
    ax_right.set_xlabel('$\mathregular{SO_2}$ oxidation rate fraction', color='k', fontsize=20, labelpad=10)
    ax_right.set_ylabel('Sulfate ($\mathregular{{\mu}g}$/$\mathregular{m^3}$)', color='k', fontsize=20, labelpad=10)
    ax_right.set_ylim(ymin, ymax)
    y_down = 0.08
    ax_right.text((xmax-xmin)*0.03+xmin, (ymax-ymin)*0.92+ymin, label_right, fontsize=22, fontweight='bold')
    ax_right.text((xmax-xmin)*0.03+xmin, (ymax-ymin)*(0.92-y_down)+ymin,  r'$\mathregular{slope:}$', fontsize=20)
    ax_right.text((xmax-xmin)*0.25+xmin, (ymax-ymin)*(0.92-y_down)+ymin,  r'$\mathregular{%.2f \pm %.2f}$' %(slope_3, pcov[0, 0]**0.5), fontsize=20)
    ax_right.text((xmax-xmin)*0.03+xmin, (ymax-ymin)*(0.85-y_down)+ymin, r'$\mathregular{intercept:}$', fontsize=20)
    ax_right.text((xmax-xmin)*0.25+xmin, (ymax-ymin)*(0.85-y_down)+ymin, r'$\mathregular{%.2f \pm %.2f}$' %(intercept_3, pcov[1, 1]**0.5), fontsize=20)
    ax_right.text((xmax-xmin)*0.03+xmin, (ymax-ymin)*(0.76-y_down)+ymin, r'$\mathregular{R^2:}$', fontsize=20)
    ax_right.text((xmax-xmin)*0.25+xmin, (ymax-ymin)*(0.76-y_down)+ymin, r'$\mathregular{%.3f}$' %(r_value_3**2), fontsize=20)
    ax_right.text((xmax-xmin)*0.03+xmin, (ymax-ymin)*(0.69-y_down)+ymin, '$\mathregular{p value:}$', fontsize=20)
    ax_right.text((xmax-xmin)*0.25+xmin, (ymax-ymin)*(0.69-y_down)+ymin, r'$\mathregular{%.4f}$' %(p_value_3), fontsize=20)

    if detail:
        print('-'*10+f'{height} m'+'-'*10)
        print(f'NOy count: {temp_dd.NOy.count()}, PAN count: {temp_dd.PAN_PPB.count()}')
        print('PAN-NOy slope:', slope_PAN, r_value_PAN**2, p_value_PAN)
        print('slope: %.1f+/-%.1f,' %(slope_3, pcov[0, 0]**0.5))
        print('intercept: %.2f+/-%.2f,' %(intercept_3, pcov[1, 1]**0.5))
        print('R2:%.3f, p value: %.3E' %(r_value_3**2, p_value_3))
        print('explicable SO4: %.2f' %(SO4_max-intercept_3))                             
        print('explicable SO4 percentage: %.1f+/-%.1f' %(explicable_SO4_percentage, pcov[1, 1]**0.5/SO4_max*100))
        print('RH: %.1f +/- %.1f' %(np.mean(RH), np.std(RH)/np.sqrt(len(RH))))
        print('SO2: %.2f +/- %.2f' %(np.mean(so2_all), np.std(so2_all)/np.sqrt(len(so2_all))))
        print('Estimated reaction time: %.1f h' %rxn_time)
        print('Max rate fraction:',np.max(slope_x[mask_3]))


# In[65]:


def plot_figure_5(aircraft_all, detail=False):
    """
    Plot Figure 5: The flight trail of data (left) utilized to calculate the correlation of
    sulfate and SO2 oxidation rate fraction (right). The altitudes of the flights are
    (a, b) 300 m, (c, d) 300 m, (e, f) 500 m, (g, h) 700 m, (i, j) 900 m, and (k, l) 1200 m, respectively. 
    """
    fig, ((ax00, ax01),
          (ax10, ax11),
          (ax20, ax21),
          (ax30, ax31),
          (ax40, ax41),
          (ax50, ax51),) = plt.subplots(6, 2, sharex=False, sharey=False, figsize=(20, 48))
    height_slope_package(aircraft_all, ax00, ax01, 20180324, 300, xmin=0, xmax=0.5, xinterval=0, box_w=20,  
                         lon_min=117.49, lon_max=131, lon_tick_num=5, 
                         lat_min=24.747, lat_max=34.5, lat_tick_num=5, 
                         ymin=2.8, ymax=5.9, lonlat_gap=1/4, label_left='(a)', label_right='(b)', detail=detail)
    height_slope_package(aircraft_all, ax10, ax11, 20180317, 300, xmin=0, xmax=0.16, xinterval=0, box_w=20,
                         lon_min=121.352, lon_max=133, lon_tick_num=4, 
                         lat_min=24.7, lat_max=33, lat_tick_num=5, 
                         ymin=0, ymax=4.0, lonlat_gap=1/4, label_left='(c)', label_right='(d)', detail=detail)
    height_slope_package(aircraft_all, ax20, ax21, 20180330, 500, xmin=0, xmax=0.12, xinterval=0, box_w=20, 
                         lon_min=132.05, lon_max=140, lon_tick_num=5, 
                         lat_min=30.15, lat_max=36, lat_tick_num=5, 
                         ymin=0, ymax=1.2, lonlat_gap=1/4, label_left='(e)', label_right='(f)', detail=detail)
    height_slope_package(aircraft_all, ax30, ax31, 20180330, 700, xmin=0, xmax=0.1, xinterval=0, box_w=6, 
                         lon_min=118.0, lon_max=124.0, lon_tick_num=5,
                         lat_min=21.775, lat_max=26, lat_tick_num=5,
                         ymin=0.1, ymax=1.9, lonlat_gap=1/4, label_left='(g)', label_right='(h)', detail=detail)
    height_slope_package(aircraft_all, ax40, ax41, 20180403, 900, xmin=0, xmax=0.04, xinterval=0, box_w=6, 
                         lon_min=118.12, lon_max=124.0, lon_tick_num=5,
                         lat_min=21.79, lat_max=26, lat_tick_num=5,
                         ymin=0.0, ymax=4.0, lonlat_gap=1/4, label_left='(i)', label_right='(j)', detail=detail)
    height_slope_package(aircraft_all, ax50, ax51, 20180403, 1200, xmin=0, xmax=0.06, xinterval=0, box_w=6, 
                         lon_min=118.04, lon_max=124.0, lon_tick_num=5,
                         lat_min=21.925, lat_max=26, lat_tick_num=5,
                         ymin=0.9, ymax=2.9, lonlat_gap=1/4, label_left='(k)', label_right='(l)', detail=detail)
    # save figure
    fig.savefig('figures/Fig.5-sulfate_SO4_oxi_rate_fraction_air_diff_area.jpg', dpi=300, bbox_inches='tight')
    plt.show()


# In[66]:


def plot_figure_6(epa_dict, aircraft_data, detail=False):
    """
    Plot Figure 6: The relationships of (a) NOx to CO at the surface and NOy to CO in the air,
    (b) SO2 to CO, and (c) SO2 to NOx at the surface and SO2 to NOy in the air. 
    """
    fig, (ax0, ax1, ax2) = plt.subplots(1, 3, sharex=False, sharey=False, figsize=(27, 6))
    temp_air = aircraft_data[aircraft_data['CO']<1000]

    label_fontsize = 24
    tick_size = 20
    index_size= 30
    label_padding = 6
    #---- CO - NOx ------
    mask = ~np.isnan(epa_dict['CO']) & ~np.isnan(epa_dict['NOx'])
    slope_epa, intercept_epa, r_value_epa, p_value_epa, std_err_epa = stats.linregress(epa_dict['CO'][mask], (epa_dict['NOx'])[mask])
    mask2 = ~np.isnan(temp_air['CO']) & ~np.isnan(temp_air['NOy'])
    slope_air, intercept_air, r_value_air, p_value_air, std_err_air = stats.linregress(temp_air['CO'][mask2], (temp_air['NOy'])[mask2])
    if detail:
        print('NOx-CO')
        print('Sfc:', slope_epa, intercept_epa, r_value_epa**2, p_value_epa, std_err_epa)
        print('Air:', slope_air, intercept_air, r_value_air**2, p_value_air, std_err_air)
    ax0.scatter(epa_dict['CO'], epa_dict['NOx'], s=25, color='deepskyblue')
    ax0.scatter(temp_air['CO'], temp_air['NOy'], s=25, color='coral')
    ax0.set_ylim(0, 260)
    ax0.set_xlim(0, 5000)
    ax0.set_xlabel('$\mathregular{CO}$ (ppb)', color='k', fontsize=label_fontsize, labelpad=label_padding)
    ax0.set_ylabel('surface $\mathregular{NO_x}$ (ppb) \n airborne $\mathregular{NO_y}$ (ppb)', color='k', fontsize=label_fontsize, labelpad=label_padding)
    ax0.xaxis.set_tick_params(labelsize=tick_size)
    ax0.yaxis.set_tick_params(labelsize=tick_size)
    (xmin, xmax), (ymin, ymax) = ax0.get_xlim(), ax0.get_ylim()
    ax0.text(xmax*0.04, ymax*0.88, '(a)', fontsize=index_size)
    ax0.plot([xmin, xmax], np.array([xmin, xmax])*slope_epa+intercept_epa, 'b')
    ax0.plot([xmin, xmax], np.array([xmin, xmax])*slope_air+intercept_air, 'r')

    #---- CO - SO2 ------
    mask = ~np.isnan(epa_dict['CO']) & ~np.isnan(epa_dict['SO2'])
    slope_epa, intercept_epa, r_value_epa, p_value_epa, std_err_epa = stats.linregress(epa_dict['CO'][mask], (epa_dict['SO2'])[mask])
    mask2 = ~np.isnan(temp_air['CO']) & ~np.isnan(temp_air['SO2'])
    slope_air, intercept_air, r_value_air, p_value_air, std_err_air = stats.linregress(temp_air['CO'][mask2], (temp_air['SO2'])[mask2])
    if detail:
        print('SO2-CO')
        print('Sfc:', slope_epa, intercept_epa, r_value_epa**2, p_value_epa, std_err_epa)
        print('Air:', slope_air, intercept_air, r_value_air**2, p_value_air, std_err_air)
    ax1.scatter(epa_dict['CO'], epa_dict['SO2'], s=25, color='deepskyblue', label='surface')
    ax1.scatter(temp_air['CO'], temp_air['SO2'], s=25, color='coral', label='airborne')
    ax1.set_ylim(0, 30)
    ax1.set_xlim(0, 5000)
    ax1.set_xlabel('$\mathregular{CO}$ (ppb)', color='k', fontsize=label_fontsize, labelpad=label_padding)
    ax1.set_ylabel('$\mathregular{SO_2}$ (ppb)', color='k', fontsize=label_fontsize, labelpad=label_padding)
    ax1.xaxis.set_tick_params(labelsize=tick_size)
    ax1.yaxis.set_tick_params(labelsize=tick_size)
    (xmin, xmax), (ymin, ymax) = ax1.get_xlim(), ax1.get_ylim()
    ax1.text(xmax*0.04, ymax*0.88, '(b)', fontsize=index_size)
    ax1.legend(ncol=2, loc='center', bbox_to_anchor=(0.5, -0.3), fontsize=26)

    #---- SO2 - NOx ------
    mask = ~np.isnan(epa_dict['SO2']) & ~np.isnan(epa_dict['NOx'])
    slope_epa, intercept_epa, r_value_epa, p_value_epa, std_err_epa = stats.linregress(epa_dict['SO2'][mask], (epa_dict['NOx'])[mask])
    mask2 = ~np.isnan(temp_air['SO2']) & ~np.isnan(temp_air['NOy'])
    slope_air, intercept_air, r_value_air, p_value_air, std_err_air = stats.linregress(temp_air['SO2'][mask2], (temp_air['NOy'])[mask2])
    if detail:
        print('NOx-SO2')
        print('Sfc:', slope_epa, intercept_epa, r_value_epa**2, p_value_epa, std_err_epa)
        print('Air:', slope_air, intercept_air, r_value_air**2, p_value_air, std_err_air)
    ax2.scatter(epa_dict['NOx'], epa_dict['SO2'], s=25, color='deepskyblue', label='surface')
    ax2.scatter(temp_air['NOy'], temp_air['SO2'], s=25, color='coral', label='airborne')
    ax2.set_xlim(0, 260)
    ax2.set_ylim(0, 30)
    ax2.set_ylabel('$\mathregular{SO_2}$ (ppb)', color='k', fontsize=label_fontsize, labelpad=label_padding)
    ax2.set_xlabel('surface $\mathregular{NO_x}$ (ppb) \n airborne $\mathregular{NO_y}$ (ppb)', color='k', fontsize=label_fontsize, labelpad=label_padding)
    ax2.xaxis.set_tick_params(labelsize=tick_size)
    ax2.yaxis.set_tick_params(labelsize=tick_size)
    (xmin, xmax), (ymin, ymax) = ax2.get_xlim(), ax2.get_ylim()
    ax2.text(xmax*0.04, ymax*0.88, '(c)', fontsize=index_size)
    # save figure
    fig.savefig('figures/Fig.6-trace_gases_relationship.jpg', dpi=300, bbox_inches='tight')
    plt.show()


# In[67]:


def plot_figure_7(slope_data, epa_loc_dict, aircraft_data):
    """
    Plot Figure 7: Correlation coefficients (R2) before and after the O3 normalization of
    CO-SO2 (a), NO-SO2 (a), and NO2-SO2 (a) of EPA stations over Western Taiwan
    illustrated in the left and right side of circles, respectively. 
    """
    # set map region
    minLat, maxLat =  21.75, 25.55
    minLon, maxLon =  119.45, 122.45
    # plot setting
    label_fontsize = 32
    rec_alpha = 0.75
    label_arg = dict(fontsize=32, weight='bold')
    source_arg = dict(marker='X', markersize=8, linestyle='None', zorder=3)
    white_box_arg = dict(edgecolor=None, facecolor='white', fill=True, alpha=0.75, zorder=20,)
    scatter_arg_sfc = dict(cmap='rainbow', s=100, alpha=1, zorder=10, linewidths=1.5)
    box_text_arg = dict(verticalalignment='center', ha='left', color='k', fontsize=20, weight='bold', zorder=30)
    region_text_arg = dict(verticalalignment='center', ha='center', color='r', fontsize=32, weight='bold')
    Taiwan_basemap_arg = dict(projection='cyl',resolution='h',llcrnrlat=minLat,
                       urcrnrlat=maxLat, llcrnrlon=minLon, urcrnrlon=maxLon)
    
    fig = plt.figure(figsize=(18, 12))
    ax1 = fig.add_axes([0.03, 0.1, 0.43, 0.8])
    ax2 = fig.add_axes([0.47, 0.1, 0.017, 0.8])
    ax3 = fig.add_axes([0.03+0.55, 0.1, 0.43, 0.8])
    
    # surface epa
    m = Basemap(ax=ax1, **Taiwan_basemap_arg)
    m.arcgisimage(service='World_Shaded_Relief', xpixels=1000, verbose= False, alpha=0.1)
    # scatter plot for plot_var
    plot_lon, plot_lat, plot_var, plot_var_err = [], [], [], []
    plot_lon_err, plot_lat_err= [], []
    for item_num in range(len(slope_data.index)):
        try:
            if slope_data['p_value_all'][item_num] & ~np.isnan(slope_data['slope_x'][item_num]):
                plot_lon.append(float(epa_loc_dict.loc['st_x', slope_data.index[item_num]]))
                plot_lat.append(float(epa_loc_dict.loc['st_y', slope_data.index[item_num]]))
                plot_var.append(slope_data['slope_x'][item_num])
                plot_var_err.append(slope_data['slope_x_err'][item_num])
            else:
                plot_lon_err.append(float(epa_loc_dict.loc['st_x', slope_data.index[item_num]]))
                plot_lat_err.append(float(epa_loc_dict.loc['st_y', slope_data.index[item_num]]))
        except:
            None
    scatter_sfc = ax1.scatter(plot_lon, plot_lat, c=plot_var, edgecolors='k', **scatter_arg_sfc, vmin=0, vmax=0.03)
    ax1.scatter(plot_lon_err, plot_lat_err, c='darkgrey', marker='X', edgecolors='k', **scatter_arg_sfc)
    clb = fig.colorbar(scatter_sfc, cax=ax2, extend='max')
    clb_label = '$\mathregular{SO_2}$ oxidation rate fraction'
    clb.set_clim(0, 0.03)  
    clb.set_label(clb_label, labelpad=14, fontsize=22)
    clb.ax.tick_params(labelsize=18)
    
    plot_lon, plot_lat, slope_list, slope_err_list = np.array(plot_lon), np.array(plot_lat), np.array(plot_var), np.array(plot_var_err)
    for region_plot in ['region_a', 'region_b', 'region_c', 'region_d', 'region_e',]:
        (p1, p2, p3, p4) = sub_region_test(24.0, 121.0, region=region_plot, boundary_only=True)
        map_region_border(p1, p2, p3, p4, ax=ax1)
        ax1.text(np.mean([p2.lon, p3.lon])-0.15, np.mean([p2.lat, p3.lat]), region_plot[-1].capitalize(),
                 **region_text_arg)
        rec_h = 0.25
        rec_w = 0.82
        ax1.add_patch(patches.Rectangle((np.mean([p1.lon, p4.lon])+0.15, np.mean([p1.lat, p4.lat])-rec_h/2+0.01),
                                        rec_w, # rectangle width
                                        rec_h, # rectangle height
                                        **white_box_arg))
        region_test = sub_region_test(plot_lat, plot_lon, region=region_plot)
        station_num = region_test.sum()
        slope_mean = np.nanmean(slope_list[region_test])
        slope_err = np.sqrt(np.nansum(slope_err_list[region_test]**2))
        ax1.text(np.mean([p1.lon, p4.lon])+0.17, np.mean([p1.lat, p4.lat]),
                 '%.4f$\mathregular{\pm}$%.4f' %(slope_mean, slope_err), **box_text_arg)
    # set label_text
    (xmin, xmax), (ymin, ymax) = ax1.get_xlim(), ax1.get_ylim()
    ax1.text((xmax-xmin)*0.03+xmin, (ymax-ymin)*0.94+ymin, '(a)', **label_arg)

    # airborne
    for region_plot in ['region_a', 'region_b', 'region_c', 'region_d', 'region_e']:
        temp_aircraft = aircraft_data[sub_region_test(aircraft_data.lat, aircraft_data.lon, region=region_plot)]
        ax3.scatter(temp_aircraft['lon'], temp_aircraft['lat'], c='b', s=3)
        (p1, p2, p3, p4) = sub_region_test(24.0, 121.0, region=region_plot, boundary_only=True)
        map_region_border(p1, p2, p3, p4, ax=ax3)
        ax3.text(np.mean([p2.lon, p3.lon])-0.15, np.mean([p2.lat, p3.lat]), region_plot[-1].capitalize(),
                 **region_text_arg)
        rec_h = 0.25
        rec_w = 0.69
        ax3.add_patch(patches.Rectangle((np.mean([p1.lon, p4.lon])+0.15, np.mean([p1.lat, p4.lat])-rec_h/2+0.01),
                                        rec_w, # rectangle width
                                        rec_h, # rectangle height
                                        **white_box_arg))
        # COO3-SO2O3
        mask_CO = ~np.isnan(temp_aircraft['SO2_over_O3']) & ~np.isnan(temp_aircraft['CO_over_O3'])
        slope_CO, intercept_CO, r_value_CO, p_value_CO, std_err_CO = stats.linregress(temp_aircraft['SO2_over_O3'][mask_CO],
                                                                                      temp_aircraft['CO_over_O3'][mask_CO])
        # NOO3-SO2O3
        mask_NO = ~np.isnan(temp_aircraft['SO2_over_O3']) & ~np.isnan(temp_aircraft['NO_over_O3'])
        slope_NO, intercept_NO, r_value_NO, p_value_NO, std_err_NO = stats.linregress(temp_aircraft['SO2_over_O3'][mask_NO],
                                                                                      temp_aircraft['NO_over_O3'][mask_NO])     
        # NO2O3-SO2O3      
        mask_NO2 = ~np.isnan(temp_aircraft['SO2_over_O3']) & ~np.isnan(temp_aircraft['NO2_over_O3'])
        slope_NO2, intercept_NO2, r_value_NO2, p_value_NO2, std_err_NO2 = stats.linregress(temp_aircraft['SO2_over_O3'][mask_NO2],
                                                                                           temp_aircraft['NO2_over_O3'][mask_NO2])
        slope_mean = 1/(slope_CO*0.16+slope_NO*24+slope_NO2*18.7+1)  
        slope_err = slope_mean*np.sqrt((0.16*std_err_CO/np.sqrt(sum(mask_CO)))**2+                                       (24*std_err_NO/np.sqrt(sum(mask_NO)))**2+                                       (18.7*std_err_NO2/np.sqrt(sum(mask_NO2)))**2)
        ax3.text(np.mean([p1.lon, p4.lon])+0.17, np.mean([p1.lat, p4.lat]),
                 '%.3f$\mathregular{\pm}$%.3f' %(slope_mean, slope_err), **box_text_arg)
    m = Basemap(ax=ax3, **Taiwan_basemap_arg)
    m.arcgisimage(service='World_Shaded_Relief', xpixels=1000, verbose= False, alpha=0.1)
    # set label
    (xmin, xmax), (ymin, ymax) = ax3.get_xlim(), ax3.get_ylim()
    ax3.text((xmax-xmin)*0.03+xmin, (ymax-ymin)*0.94+ymin, '(b)', **label_arg)
    # save figure
    fig.savefig('figures/Fig.7-Map_region_SO2_oxi_rate_frac_sfc_air.jpg', dpi=300, bbox_inches='tight')
    plt.show()


# In[68]:


def plot_figure_S1(aircraft_all):
    """
    Plot Figure S1: Flight tracks of EMeRGe-Asia campaign over West Pacific.
    """
    aircraft_all = aircraft_all_input()
    fig = plt.figure(figsize=(12, 9))
    ax1 = fig.add_axes([0.05, 0.05, 0.95, 0.95])

    # set map region
    minLat, maxLat = 10., 40.
    minLon, maxLon = 110., 140.

    date_plot = list(sorted(set(aircraft_all.date)))[2:-1]
    c_l = [plt.cm.rainbow(i/11) for i in range(len(date_plot))]
    for i in range(len(date_plot)):
        date = date_plot[i]
        temp = aircraft_all[aircraft_all.date==date]
        ax1.plot(temp.lon, temp.lat, color=c_l[i], linewidth=1.5, label=str(int(date)))

    ax1.legend(loc='upper left', borderaxespad=1, fontsize=16, edgecolor='white', frameon=True)

    m = Basemap(projection='cyl',resolution='h',llcrnrlat=minLat,         urcrnrlat = maxLat, llcrnrlon = minLon, urcrnrlon = maxLon, ax=ax1)
    m.arcgisimage(service='World_Shaded_Relief', xpixels=1000, verbose= False, alpha=0.1)
    # set parallels and meridians label
    parallels = np.linspace(minLat, maxLat, num=7, endpoint=True)
    meridians = np.linspace(minLon, maxLon, num=7, endpoint=True)
    ax1.set_xticks(meridians)
    ax1.set_yticks(parallels)
    ax1.set_xticklabels(['%.1f$\degree$E' %i for i in meridians])
    ax1.set_yticklabels(['%.1f$\degree$N' %i for i in parallels])
    ax1.xaxis.set_tick_params(labelsize=20, pad=10)
    ax1.yaxis.set_tick_params(labelsize=20, pad=10)
    plt.show()
    # save figure
    fig.savefig('figures/Fig.S1-Map_Of_Flight_Track_over_West_Pacific.jpg', dpi=300, bbox_inches='tight')
    plt.show()


# In[69]:


def plot_figure_S2(pm_sinica, epa_loc_dict, p_to_en):
    """
    Plot Figure S2: The geographic locations of the 11 sampling sites the averaged composition of PM2.5.
    """
    fig = plt.figure(figsize=(10, 18))
    ax1 = fig.add_axes([0.1, 0.1, 0.7, 0.8])
    # plot map
    minLat, maxLat =  21.75, 25.55
    minLon, maxLon =  119.45, 122.45
    m = Basemap(projection='cyl', resolution='h',                llcrnrlat=minLat, urcrnrlat=maxLat,                llcrnrlon=minLon, urcrnrlon=maxLon, ax=ax1)
    m.arcgisimage(service='World_Shaded_Relief', xpixels = 1000, verbose= False, alpha=0.1)
    source_arg = dict(marker='X', markersize=8, linestyle='None', zorder=3)

    plot_place = sorted(set(pm_sinica.place))
    c_l = [plt.cm.rainbow(i/11) for i in range(11)]
    for i in range(11):
        ax1.scatter(float(epa_loc_dict.loc['st_x', plot_place[i]]),
                    float(epa_loc_dict.loc['st_y', plot_place[i]]),
                    color=c_l[i], s=100, edgecolors='k', label=p_to_en[plot_place[i]])
    ax1.legend(loc='center left', bbox_to_anchor=(0.99, 0.5), fontsize=20)
    # save figure
    fig.savefig('figures/Fig.S2-station_Map.jpg', dpi=300, bbox_inches='tight')
    plt.show()


# In[70]:


def plot_figure_S3(aircraft_data, k_list):
    """
    Plot Figure S3: elationship between VOCs-OH radical reaction rates and the summation of
    other trace gases (CH4, CO, NO, NO2, and SO2)-OH reaction rates. Reaction rate constant “k” of
    each trace gas is the second-order or pseudo-second-order reaction rate constant (kA[A][OH],
    which A represents trace gases we study in this research) at 298K. VOCs include formaldehyde, methanol,
    acetonitrile, acetaldehydem acetone, isoprene, benzene, toluene, and xylene.
    """
    # plot setting parameters
    label_fontsize = 22
    legend_dict = dict(ncol=2, bbox_to_anchor=(0.5, -0.22), loc='center', borderaxespad=0, labelspacing=0.5, frameon=True, fontsize=28)
    epa_scatter_s = 100
    epa_scatter_alpha = 0.5
    ytick_size=18
    textsize = 18
    
    data_select_1 = np_or(sub_region_test(np.array(aircraft_data['lat']), np.array(aircraft_data['lon']), region='region_a'), sub_region_test(np.array(aircraft_data['lat']), np.array(aircraft_data['lon']), region='region_b'))
    data_select_2 = np_or(sub_region_test(np.array(aircraft_data['lat']), np.array(aircraft_data['lon']), region='region_c'), sub_region_test(np.array(aircraft_data['lat']), np.array(aircraft_data['lon']), region='region_d'))
    data_select_3 = np_or(data_select_1, data_select_2)
    data_select_4 = np_or(sub_region_test(np.array(aircraft_data['lat']), np.array(aircraft_data['lon']), region='region_e'), data_select_3)
    temp_aircraft = aircraft_data[data_select_4 == True]

    x_list = temp_aircraft['CH4']*1000*k_list['CH4'] + temp_aircraft['CO']*k_list['CO'] +             temp_aircraft['NO']*k_list['NO'] + temp_aircraft['NO2']*k_list['NO2'] + temp_aircraft['SO2']*k_list['SO2']
    y_list = np.zeros(len(temp_aircraft))
    VOC_list = ['FOR', 'MET', 'ACN', 'ACA', 'ACE', 'ISO', 'BEN', 'TOL', 'XYL']
    for i in range(len(VOC_list)):
        y_list += temp_aircraft[VOC_list[i]]/1000*k_list[VOC_list[i]]
    x_input = x_list*2.46e10
    y_input = y_list*2.46e10

    label_CO = '$\mathregular{k_{CO}}$[$\mathregular{CO}$]'
    label_NO = '$\mathregular{k_{NO}}$[$\mathregular{NO}$]'
    label_NO2 = '$\mathregular{k_{NO_2}}$[$\mathregular{NO_2}$]'
    label_SO2 = '$\mathregular{k_{SO_2}}$[$\mathregular{SO_2}$]'
    label_CH4 = '$\mathregular{k_{CH_4}}$[$\mathregular{CH_4}$]'
    x_label = '{}+{}+{}+{}+{}'.format(label_CH4, label_CO, label_NO, label_NO2, label_SO2)
    y_label = '$\mathregular{k_{VOCs}}$[$\mathregular{VOCs}$]'
    
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
    ax.scatter(x_input, y_input, color='k')

    mask = ~np.isnan(x_input) & ~np.isnan(y_input)
    slope, intercept, r_value, p_value, std_err = stats.linregress(x_input[mask], y_input[mask])       

    xmin, xmax = 0, x_input.max()*1.1
    ymin, ymax = y_input.min()/1.1, y_input.max()*1.1
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    # plot regression line
    ax.plot([xmin, xmax], [xmin*slope+intercept, xmax*slope+intercept], 'r')
    # label and text setting
    ax.set_xlabel('%s ($\mathregular{s^{-1}}$)'%(x_label), color='k', fontsize = label_fontsize, labelpad=10)
    ax.set_ylabel('%s ($\mathregular{s^{-1}}$)'%(y_label), color='k', fontsize = label_fontsize, labelpad=10)
    ax.xaxis.set_tick_params(labelsize=ytick_size)
    ax.yaxis.set_tick_params(labelsize=ytick_size)
    ax.text((xmax-xmin)*0.65+xmin, (ymax-ymin)*0.23+ymin, 'y = %.3fx + %.3f' %(slope, intercept), fontsize=18)
    ax.text((xmax-xmin)*0.65+xmin, (ymax-ymin)*0.15+ymin, '$\mathregular{R^2}$= %.3f' %(r_value**2), fontsize=18)
    p_value_print = r'p-value=$\mathregular{%.3f}$' %(p_value) if p_value > 0.001 else r'p-value=$\mathregular{%.1e}$' %(p_value)
    ax.text((xmax-xmin)*0.65+xmin, (ymax-ymin)*0.07+ymin, p_value_print, fontsize=textsize)
    # save figure
    fig.savefig('figures/Fig.S3-VOC_rate-others_rate.jpg', dpi=300, bbox_inches='tight')
    plt.show()


# In[71]:


def plot_figure_S4(slope_data, epa_loc_dict, detail=False):
    """
    Plot Figure S4: Correlation coefficients (R2) before and after the O3 normalization of
    CO-SO2 (a), NO-SO2 (a), and NO2-SO2 (a) of EPA stations over Western Taiwan
    illustrated in the left and right side of circles, respectively. 
    """
    fig = plt.figure(figsize=(27, 16))
    ax1 = fig.add_axes([0.03, 0.1, 0.28, 0.8])
    ax2 = fig.add_axes([0.03+0.3, 0.1, 0.28, 0.8])
    ax3 = fig.add_axes([0.03+0.6, 0.1, 0.28, 0.8])
    ax4 = fig.add_axes([0.935, 0.11, 0.017, 0.78])

    ax_list = (ax1, ax2, ax3)
    var_list = ['CO', 'NO', 'NO2']
    label_text = ['a', 'b', 'c']
    
    # map region
    min_lat, max_lat = 22.25, 24.75            
    min_lon, max_lon = 119.5, 121
    # plot setting
    scatter_arg_sfc = dict(markersize=20, alpha=1, zorder=6, linestyle='None', marker='o', fillstyle='right')
    label_arg = dict(fontsize=32, weight='bold')
    label_fontsize = 32

    for index in range(3):
        """plot variable and wind on the map"""
        col = var_list[index]
        ax1 = ax_list[index]
        # plot map
        m = Basemap(projection='cyl', resolution='h',                    llcrnrlat=min_lat, urcrnrlat=max_lat,                    llcrnrlon=min_lon, urcrnrlon=max_lon, ax=ax1)
        # select Basemap background
        m.arcgisimage(service='World_Shaded_Relief', xpixels=1000, verbose= False, alpha=0.1)
        # scatter plot
        plot_lon, plot_lat, plot_var, plot_var_2, in_map_region = [], [], [], [], []
        for item_num in range(len(slope_data)):
            lon_item = float(epa_loc_dict.loc['st_x', slope_data.index[item_num]])
            lat_item = float(epa_loc_dict.loc['st_y', slope_data.index[item_num]])
            lon_in_map = lon_item>=min_lon and lon_item<=max_lon
            lat_in_map = lat_item>=min_lat and lat_item<=max_lat
            in_map_region.append(lon_in_map and lat_in_map)
            try:
                plot_lon.append(float(epa_loc_dict.loc['st_x', slope_data.index[item_num]]))
                plot_lat.append(float(epa_loc_dict.loc['st_y', slope_data.index[item_num]]))
                plot_var.append(slope_data[f'{col}/O3-SO2/O3_R2'][item_num])
                plot_var_2.append(slope_data[f'{col}-SO2_R2'][item_num])
            except:
                None
        for item_num in range(len(plot_var)):
            ax1.plot(plot_lon[item_num], plot_lat[item_num], c=plt.cm.rainbow(plot_var[item_num]),
                     markerfacecoloralt=plt.cm.rainbow(plot_var_2[item_num]), markeredgecolor='k', **scatter_arg_sfc)              
        # set label_text
        xmin, xmax = ax1.get_xlim()
        ymin, ymax = ax1.get_ylim()
        ax1.text((xmax-xmin)*0.03+xmin, (ymax-ymin)*0.95+ymin, f'({label_text[index]})', **label_arg)
        if detail:
            print(f'--------   {col}   --------')
            R2_before_O3_mean = slope_data[f'{col}-SO2_R2'][in_map_region].mean()
            R2_before_O3_err = slope_data[f'{col}-SO2_R2'][in_map_region].std()/np.sqrt(slope_data[f'{col}-SO2_R2'][in_map_region].count())
            R2_after_O3_mean = slope_data[f'{col}/O3-SO2/O3_R2'][in_map_region].mean()
            R2_after_O3_err = slope_data[f'{col}/O3-SO2/O3_R2'][in_map_region].std()/np.sqrt(slope_data[f'{col}/O3-SO2/O3_R2'][in_map_region].count())
            print(f'slope R2 before O3 mean: {R2_before_O3_mean:.3f}+/-{R2_before_O3_err:.3f}')
            print(f'slope R2 before O3 mean: {R2_after_O3_mean:.3f}+/-{R2_after_O3_err:.3f}')

    # set colorbar
    scatter_sfc = ax1.scatter(plot_lon, plot_lat, c=plot_var, edgecolors='k', cmap='rainbow', s=0, alpha=1, vmin=0.2, vmax=0.8)
    norm = matplotlib.colors.Normalize(vmin=0.2, vmax=0.8)
    clb = fig.colorbar(scatter_sfc, cax=ax4, extend='both')
    clb_label = '$\mathregular{R^2}$'
    clb.set_clim(0.2, 0.8)
    clb.set_label(clb_label, labelpad=14, fontsize=28)
    clb.ax.tick_params(labelsize=22)
    # save figure
    fig.savefig('figures/Fig.S4-Map_R2_both_slope_after.jpg', dpi=300, bbox_inches='tight')
    plt.show()


# In[72]:


def plot_figure_S5(aircraft_data):
    """
    Plot Figure S5: Five self-defined regions over Western Taiwan with the flight tracks in these regions.
    """

    # set map region
    minLat, maxLat =  21.75, 25.55
    minLon, maxLon =  119.45, 122.45
    Taiwan_basemap_arg = dict(projection='cyl',resolution='h',llcrnrlat=minLat,
                       urcrnrlat=maxLat, llcrnrlon=minLon, urcrnrlon=maxLon)
    
    fig = plt.figure(figsize=(14, 10.5))
    ax1 = fig.add_axes([0.05, 0.05, 0.95, 0.95])
    for region_plot in ['region_a', 'region_b', 'region_c', 'region_d', 'region_e']:
        temp_aircraft = aircraft_data[sub_region_test(aircraft_data.lat, aircraft_data.lon, region=region_plot)]
        ax1.scatter(temp_aircraft['lon'], temp_aircraft['lat'], c='b', s=3)
        (p1, p2, p3, p4) = sub_region_test(24.0, 121.0, region=region_plot, boundary_only=True)
        map_region_border(p1, p2, p3, p4, ax=ax1)
        ax1.text(np.mean([p2.lon, p3.lon])-0.15, np.mean([p2.lat, p3.lat]), region_plot[-1].capitalize(),
                 verticalalignment='center', ha='center', color='r', fontsize=32, weight='bold')
    m = Basemap(ax=ax1, **Taiwan_basemap_arg)
    m.arcgisimage(service='World_Shaded_Relief', xpixels=1000, verbose= False, alpha=0.1)
    # save figure
    fig.savefig('figures/Fig.S5-Map_R2_both_slope_after.jpg', dpi=300, bbox_inches='tight')
    plt.show()


# In[73]:


def plot_figure_S6(slope_data, epa_dict, epa_loc_dict, aircraft_data):
    """
    Plot Figure S6: The average wind speed of each region at the surface (a) and in the air (b). 
    """
    fig = plt.figure(figsize=(18, 12))
    ax1 = fig.add_axes([0.03, 0.1, 0.43, 0.8])
    ax2 = fig.add_axes([0.47, 0.1, 0.017, 0.8])
    ax3 = fig.add_axes([0.03+0.55, 0.1, 0.43, 0.8])
    
    # set map region
    minLat, maxLat =  21.75, 25.55
    minLon, maxLon =  119.45, 122.45
    # plot setting
    scatter_arg_sfc = dict(cmap='rainbow', s=100, alpha=1, zorder=10, linewidths=1.5)
    label_arg = dict(fontsize=32, weight='bold')
    source_arg = dict(marker='X', markersize=8, linestyle='None', zorder=3)
    white_box_arg = dict(edgecolor=None, facecolor='white', fill=True, alpha=0.75, zorder=20,)
    box_text_arg = dict(verticalalignment='center', ha='left', color='k', fontsize=20, weight='bold', zorder=30)
    region_text_arg = dict(verticalalignment='center', ha='center', color='r', fontsize=32, weight='bold')
    Taiwan_basemap_arg = dict(projection='cyl',resolution='h',llcrnrlat=minLat,
                       urcrnrlat=maxLat, llcrnrlon=minLon, urcrnrlon=maxLon)
    
    # surface epa
    m = Basemap(ax=ax1, **Taiwan_basemap_arg)
    m.arcgisimage(service='World_Shaded_Relief', xpixels=1000, verbose= False, alpha=0.1)
    # scatter plot for plot_var
    plot_lon, plot_lat, plot_var, plot_var_err = [], [], [], []
    plot_lon_err, plot_lat_err= [], []
    u_list, v_list = [], []
    for item_num in range(len(slope_data.index)):
        try:
            if slope_data['p_value_all'][item_num] & ~np.isnan(slope_data['slope_x'][item_num]):
                plot_lon.append(float(epa_loc_dict.loc['st_x', slope_data.index[item_num]]))
                plot_lat.append(float(epa_loc_dict.loc['st_y', slope_data.index[item_num]]))
                plot_var.append(slope_data['slope_x'][item_num])
                plot_var_err.append(slope_data['slope_x_err'][item_num])
                temp_epa = epa_dict[epa_dict.station == slope_data.index[item_num]]
                u, v = wind_calc.wind_spddir_to_uv(temp_epa['WIND_SPEED'], temp_epa['WIND_DIREC'])
                u_list.append(np.nanmean(u))
                v_list.append(np.nanmean(v))
            else:
                plot_lon_err.append(float(epa_loc_dict.loc['st_x', slope_data.index[item_num]]))
                plot_lat_err.append(float(epa_loc_dict.loc['st_y', slope_data.index[item_num]]))
        except:
            None
    scatter_sfc = ax1.scatter(plot_lon, plot_lat, c=plot_var, edgecolors='k', **scatter_arg_sfc)
    ax1.scatter(plot_lon_err, plot_lat_err, c='darkgrey', marker='X', edgecolors='k', vmin=0, vmax=0.03, **scatter_arg_sfc)
    clb = fig.colorbar(scatter_sfc, cax=ax2, extend='max')
    clb_label = '$\mathregular{SO_2}$ oxidation rate fraction'
    clb.set_clim(0, 0.03)  
    clb.set_label(clb_label, labelpad=14, fontsize=22)
    clb.ax.tick_params(labelsize=18)
    
    plot_lon, plot_lat, u_list, v_list = np.array(plot_lon), np.array(plot_lat), np.array(u_list), np.array(v_list)
    for region_plot in ['region_a', 'region_b', 'region_c', 'region_d', 'region_e',]:
        (p1, p2, p3, p4) = sub_region_test(24.0, 121.0, region=region_plot, boundary_only=True)
        map_region_border(p1, p2, p3, p4, ax=ax1)
        ax1.text(np.mean([p2.lon, p3.lon])-0.15, np.mean([p2.lat, p3.lat]), region_plot[-1].capitalize(),
                 **region_text_arg)
        rec_h = 0.3
        rec_w = 0.6
        ax1.add_patch(patches.Rectangle((np.mean([p1.lon, p4.lon])+0.15, np.mean([p1.lat, p4.lat])-rec_h/2+0.01),
                                        rec_w, # rectangle width
                                        rec_h, # rectangle height
                                        **white_box_arg))
        region_test = sub_region_test(plot_lat, plot_lon, region=region_plot)
        station_num = region_test.sum()
        u_mean = u_list[region_test].mean()
        u_err = u_list[region_test].std()/np.sqrt(station_num)
        v_mean = v_list[region_test].mean()
        v_err = v_list[region_test].std()/np.sqrt(station_num)
        ax1.text(np.mean([p1.lon, p4.lon])+0.17, np.mean([p1.lat, p4.lat])+rec_h/4,
                 'u: %.1f$\mathregular{\pm}$%.1f' %(u_mean, u_err), **box_text_arg)
        ax1.text(np.mean([p1.lon, p4.lon])+0.17, np.mean([p1.lat, p4.lat])-rec_h/4,
                 'v: %.1f$\mathregular{\pm}$%.1f' %(v_mean, v_err), **box_text_arg)
    # set label_text
    (xmin, xmax), (ymin, ymax) = ax1.get_xlim(), ax1.get_ylim()
    ax1.text((xmax-xmin)*0.03+xmin, (ymax-ymin)*0.94+ymin, '(a)', **label_arg)

    # airborne
    for region_plot in ['region_a', 'region_b', 'region_c', 'region_d', 'region_e']:
        temp_aircraft = aircraft_data[sub_region_test(aircraft_data.lat, aircraft_data.lon, region=region_plot)]
        ax3.scatter(temp_aircraft['lon'], temp_aircraft['lat'], c='b', s=3)
        (p1, p2, p3, p4) = sub_region_test(24.0, 121.0, region=region_plot, boundary_only=True)
        map_region_border(p1, p2, p3, p4, ax=ax3)
        ax3.text(np.mean([p2.lon, p3.lon])-0.15, np.mean([p2.lat, p3.lat]), region_plot[-1].capitalize(),
                 **region_text_arg)
        rec_h = 0.3
        rec_w = 0.6
        ax3.add_patch(patches.Rectangle((np.mean([p1.lon, p4.lon])+0.15, np.mean([p1.lat, p4.lat])-rec_h/2+0.01),
                                        rec_w, # rectangle width
                                        rec_h, # rectangle height
                                        **white_box_arg))
        u_mean = temp_aircraft['u_ma'].mean()
        u_err = temp_aircraft['u_ma'].std()/np.sqrt(len(temp_aircraft))
        v_mean = temp_aircraft['v_ma'].mean()
        v_err = temp_aircraft['v_ma'].std()/np.sqrt(len(temp_aircraft))
        ax3.text(np.mean([p1.lon, p4.lon])+0.17, np.mean([p1.lat, p4.lat])+rec_h/4,
                'u: %.1f$\mathregular{\pm}$%.1f' %(u_mean, u_err), **box_text_arg)
        ax3.text(np.mean([p1.lon, p4.lon])+0.17, np.mean([p1.lat, p4.lat])-rec_h/4,
                'v: %.1f$\mathregular{\pm}$%.1f' %(v_mean, v_err), **box_text_arg)
    m = Basemap(ax=ax3, **Taiwan_basemap_arg)
    m.arcgisimage(service='World_Shaded_Relief', xpixels=1000, verbose= False, alpha=0.1)
    # set label
    (xmin, xmax), (ymin, ymax) = ax3.get_xlim(), ax3.get_ylim()
    ax3.text((xmax-xmin)*0.03+xmin, (ymax-ymin)*0.94+ymin, '(b)', **label_arg)
    # save figure
    fig.savefig('figures/Fig.S6-Map_region_wind_speed_sfc_air.jpg', dpi=300, bbox_inches='tight')
    plt.show()


# In[76]:


def main():
    # check if figures directory exists
    if not os.path.exists('figures'):
        cmd = 'mkdir figures'
        os.popen(cmd)
    # dataset input
    epa_dict, March_epa_dict, epa_loc_dict = epa_data_input()
    slope_data = slope_data_propcess(epa_dict)
    pm_sinica = pm_sinica_input()
    aircraft_data = aircraft_data_input()
    aircraft_all = aircraft_all_input()
    # constant dictionary
    p_to_en = {'二林':'Erlin', '大里':'Dali', '沙鹿':'Shalu', '麥寮':'Mailiao',
               '豐原':'Fengyuan', '線西':'Xianxi', '南投':'Nantou', '竹山':'Chushan', '埔里':'Puli',
               '忠明':'Chungming', '斗六':'Douliou'}
    k_list = {'FOR':9.37e-12, 'MET':9.44e-13, 'ACN':3.00e-14, 
              'ACA':1.58e-11, 'ACE':2.19e-13, 'ISO':1.01e-10, 
              'BEN':1.22e-12, 'TOL':5.63e-12, 'XYL':2.30e-11, 
              'DMS':5.00e-12, 'MEK':1.15e-12,
              'CO':2.40e-13,  'NO':3.60e-11,  'NO2':2.80e-11, 
              'CH4':6.30e-15, 'SO2':1.5e-12} 
    # plot paper figures
    print(plot_figure_2.__doc__)
    plot_figure_2(epa_dict)
    print(plot_figure_3.__doc__)
    plot_figure_3(pm_sinica, March_epa_dict, p_to_en, aircraft_data, detail=False)
    print(plot_figure_4.__doc__)
    plot_figure_4(pm_sinica, March_epa_dict, p_to_en)
    print(plot_figure_5.__doc__)
    plot_figure_5(aircraft_all, detail=False)
    print(plot_figure_6.__doc__)
    plot_figure_6(epa_dict, aircraft_data, detail=False)
    print(plot_figure_7.__doc__)
    plot_figure_7(slope_data, epa_loc_dict, aircraft_data)
    print(plot_figure_S1.__doc__)
    plot_figure_S1(aircraft_all)
    print(plot_figure_S2.__doc__)
    plot_figure_S2(pm_sinica, epa_loc_dict, p_to_en)
    print(plot_figure_S3.__doc__)
    plot_figure_S3(aircraft_data, k_list)
    print(plot_figure_S4.__doc__)
    plot_figure_S4(slope_data, epa_loc_dict, detail=False)
    print(plot_figure_S5.__doc__)
    plot_figure_S5(aircraft_data)
    print(plot_figure_S6.__doc__)
    plot_figure_S6(slope_data, epa_dict, epa_loc_dict, aircraft_data)


# In[77]:


if __name__ == '__main__':
    main()


# In[ ]:





# In[ ]:




