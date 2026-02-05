"""
Script to produce interpolations used to compare ARCADE year-long results with
BGR/PMCC data.

by Rodrigo De Negri
rodrigo.denegri@uca.fr
2025-01-24

platform: Ubuntu 24.04
"""

from os.path import join
import sys
import pickle
from pandas import read_csv
import numpy as np
import pandas as pd
from scipy import interpolate


def interp_this(doy_lst, azdev_lst, smooth):
    """Function to interpolate using doys and azdevs"""
    to_remove = []
    for i, azdev_i in enumerate(azdev_lst):
        # Check for nans
        if np.isnan(azdev_i):
            to_remove.append(i)
        # Check for empty values
        if not azdev_i and azdev_i != 0:
            to_remove.append(i)
    doy_lst_cp = []
    azdev_lst_cp = []
    print(f"[interp_this] Note: I'm removing {len(to_remove)} Nan values...")
    j = 0
    for i in to_remove:
        for l in doy_lst[j:i]:
            doy_lst_cp.append(l)
        for l in azdev_lst[j:i]:
            azdev_lst_cp.append(l)
        j = i+1
    if j < 366:
        for l in doy_lst[j:]:
            doy_lst_cp.append(l)
        for l in azdev_lst[j:]:
            azdev_lst_cp.append(l)

    interp_vals = interpolate.splrep(doy_lst_cp, azdev_lst_cp, s=smooth)
    min_doy = int(doy_lst_cp[0])
    max_doy = int(doy_lst_cp[-1])
    new_doys = np.linspace(min_doy, max_doy, 3*(max_doy-min_doy))
    new_azdevs = interpolate.splev(new_doys, interp_vals, der=0)
    return new_doys, new_azdevs

def interpolate_case(data, case, stat, smooth, data_dump):
    """
    Function to create the list wihth the interpolation related results
    from the back-azimuth table
    """
    doy_all = []
    doy_all_std = []
    av_bazdev_all = []
    av_bazdev_all_std = []
    av_bazdev_all_std_only_4 = []
    bazdev_s = []
    n_times = []
    print(f"[interpolate_case] Station {stat}")
    cond_sta = data['StaNam'] == stat
    data = data[cond_sta]
    for doy in range(1, 366):
        data_doy = data[data['DOY'] == doy]
        num_times = len(data_doy)
        n_times.append(num_times)
        if num_times == 0:
            av_bazdev_all.append(np.nan)
        else:
            baz_a = data_doy['BazDevA']
            baz_t = data_doy['BazDevT']
            baz_s = data_doy['BazDevS']
            if np.isnan(baz_a.mean()):
                if np.isnan(baz_t.mean()):
                    if np.isnan(baz_s.mean()):
                        print(baz_s)
                        breakpoint()
                    else:
                        av_bazdev_all.append(baz_s.mean())
                        av_bazdev_all_std.append(baz_s.std())
                else:
                    av_bazdev_all.append(baz_t.mean())
                    av_bazdev_all_std.append(baz_t.std())
            else:
                av_bazdev_all.append(baz_a.mean())
                av_bazdev_all_std.append(baz_a.std())

        # count how many of the four times are present
        if len(data_doy) == 4:
            doy_all_std.append(doy)
            av_bazdev_all_std_only_4.append(data_doy['BazDevA'].std())

        strato_arriv_exists = False
        for bazdev_s_i in data_doy['#BazDS']:
            if bazdev_s_i > 0:
                strato_arriv_exists = True
        doy_all.append(doy)
        bazdev_s.append(strato_arriv_exists)

    doy_interp, av_bazdev_interp = interp_this(
        doy_all,
        av_bazdev_all,
        smooth
    )

    # interpolate the standard deviations
    doy_interp_std = []
    av_std_interp = []
    if len(av_bazdev_all_std_only_4) > 10:
        doy_interp_std, av_std_interp = interp_this(
                doy_all_std,
                av_bazdev_all_std_only_4,
                smooth
            )

    # Export data in pickle
    # NOTE: this is an ad-hoc dictionary (similar to what is in
    #       get_stat_name function in plot_example.py)
    stat_dic = {
        'IS02': 'I02AR',
        'IS08': 'I08BO',
        'IS41': 'I41PY',
        }
    pickle_out = f"{case}_{stat_dic[stat]}.pkl"
    pickle.dump(
        [
            doy_interp,
            av_bazdev_interp,
            np.asarray(doy_all),
            np.asarray(av_bazdev_all),
            np.asarray(av_bazdev_all_std),
            bazdev_s,
            np.asarray(n_times),
            doy_interp_std,
            av_std_interp,
            smooth,
        ],
        open(join(data_dump, pickle_out), 'wb')
    )
    print(f"-> SAVED: {join(data_dump, pickle_out)}")

if __name__ == '__main__':
    print("\n=====================")
    print("Running interpolation")
    print("=====================")

    # ========
    # Settings
    # ========
    header_names = {'Year':int,
                    'DOY':int, 
                    'Seconds':int,
                    'SouNum':int, 
                    'StaNum':int, 
                    'StaNam':str,
                    'TrueBaz':float, 
                    'BazDevS':float,
                    '#BazDS':float,
                    'BazDevT':float, 
                    '#BazDT':float , 
                    'BazDevA':float,
                    'StdBDA':float, 
                    'Ill':bool}

    # Data
    INTERP_OUTPUT_PATH = "."
    print(f"Location of output: <{INTERP_OUTPUT_PATH}>")

    # Input call for data loading and interpolating
    CASE = 'Puye'
    SMOOTH = 10
    print(f"Case: <{CASE}>")
    print(f"Smoothness: <{SMOOTH}>")
    
    # Predefined project names, output of ARCADE
    # Note: all these runs are strato & thermo, adaptive perturbations
    #       of climatologies
    tab_CC = read_csv('azimuth_deviation_table.txt',
                           sep='\\s+', header=0,
                           names=[i for i in header_names])
    tab_CC = tab_CC.drop_duplicates()
    tab_CC = tab_CC[tab_CC['Ill'] == False]
    # ===============================================================
    # TEMPORARY FIX, ADDRESS THIS ISSUE IN ARCADE
    # Flip back-azimuth deviation if value is higher than 10 degrees
    #================================================================
    cond0 = tab_CC['StaNam'] == 'IS02'
    cond2 = tab_CC['DOY'] < 250
    cond1 = tab_CC['DOY'] > 150

    table_filt = tab_CC[cond0 & cond1 & cond2]

    for baz in ['BazDevS', 'BazDevT', 'BazDevA']:
        table_filt_aux = table_filt[table_filt[baz]>14]
        for ind in table_filt_aux.index:
            tab_CC.at[ind, baz] = -table_filt_aux.loc[ind][baz]
    #////////////////////////////////////////////////////////////////

    interpolate_case(tab_CC, CASE, 'IS02', SMOOTH, INTERP_OUTPUT_PATH)
    interpolate_case(tab_CC, CASE, 'IS08', SMOOTH, INTERP_OUTPUT_PATH)
    interpolate_case(tab_CC, CASE, 'IS41', SMOOTH, INTERP_OUTPUT_PATH)
