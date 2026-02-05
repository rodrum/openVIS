"""
Example script to plot VIS results.

by Rodrigo De Negri
rodrigo.denegri@uca.fr
2025-01-24

platform: Ubuntu 24.04
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys
from datetime import datetime, timedelta, timezone
import toml
import os
from os.path import join
from obspy.geodetics.base import gps2dist_azimuth
import matplotlib.dates as mdates
import matplotlib.patches as mpatches

def get_stat_name(stat_name):
    """
    Temporary patch to find interpolations following database names of stations
    Used for back-azimuth deviation based tolerances.
    """
    stat_dic = {
        'IS01': 'I01AR',
        'IS02': 'I02AR',
        'IS03': 'I03AU',
        'IS04': 'I04AU',
        'IS05': 'I05AU',
        'IS06': 'I06AU',
        'IS07': 'I07AU',
        'IS08': 'I08BO',
        'IS09': 'I09BR',
        'IS10': 'I10CA',
        'IS11': 'I11CV',
        'IS12': 'I12CR',
        'IS13': 'I13CL',
        'IS14': 'I14CL',
        'IS15': 'I15CN',
        'IS16': 'I16CN',
        'IS17': 'I17CI',
        'IS18': 'I18DK',
        'IS19': 'I19DJ',
        'IS20': 'I20EQ',
        'IS21': 'I21FR',
        'IS22': 'I22FR',
        'IS23': 'I23FR',
        'IS24': 'I24FR',
        'IS25': 'I25FR',
        'IS26': 'I26DE',
        'IS27': 'I27DE',
        'IS29': 'I29IR',
        'IS30': 'I30JP',
        'IS31': 'I31KZ',
        'IS32': 'I32KE',
        'IS33': 'I33MG',
        'IS34': 'I34MN',
        'IS35': 'I35NA',
        'IS36': 'I36NZ',
        'IS37': 'I37NO',
        'IS38': 'I38PA',
        'IS39': 'I39PW',
        'IS40': 'I40PN',
        'IS41': 'I41PY',
        'IS42': 'I42PT',
        'IS43': 'I43RU',
        'IS44': 'I44RU',
        'IS45': 'I45RU',
        'IS46': 'I46RU',
        'IS47': 'I47ZA',
        'IS48': 'I48TN',
        'IS49': 'I49GB',
        'IS50': 'I50GB',
        'IS51': 'I51GB',
        'IS52': 'I52GB',
        'IS53': 'I53US',
        'IS54': 'I54US',
        'IS55': 'I55US',
        'IS56': 'I56US',
        'IS57': 'I57US',
        'IS58': 'I58US',
        'IS59': 'I59US',
        'IS60': 'I60US',
    }
    if stat_name in stat_dic:
        return stat_dic[stat_name]
    else:
        return stat_name


def get_volcanoes_from_gvp_database(volcano_path: str) -> pd.DataFrame:
    """
    Get a list of interesting volcanoes from the very complete GVP database.

    Only volcanoes that have a recorded eruption since 2003 are selected in the
    VIS database.

    Returns
    -------
    :class:`pandas.DataFrame`
        Filtered dataframe containing the volcanoes.

    Notes
    -----
    The GVP database can be downloaded at:
    https://volcano.si.edu/volcanolist_holocene.cfm.
    """

    volcanoes_all_data = pd.read_csv(
        volcano_path,
        sep=";",
        encoding="latin-1",
        dtype={'Volcano Number': int,
               'Volcano Name': str,
               'Country': str,
               'Last Known Eruption': str,
               'Latitude': float,
               'Longitude': float,
               'Elevation (m)': float,
               },
        decimal=','
    )

    volcanoes = volcanoes_all_data[
       [
            "Volcano Number",
            "Volcano Name",
            "Country",
            "Last Known Eruption",
            "Latitude",
            "Longitude",
            "Elevation (m)",
        ]
    ].copy()
    lke = "Last Known Eruption"
    # Filter all volcanoes with lke unknown and BCE (BC equivalent)
    volcanoes.drop(
        volcanoes[
            (volcanoes[lke] == "Unknown") | (volcanoes[lke].str.contains("BCE"))
        ].index,
        inplace=True,
    )
    # Filter all volcanoes with lke before 2003
    # (year where we have most of the IMS data)
    volcanoes.drop(
        volcanoes[
            volcanoes[lke].str.extract(r"([\d]+)", expand=False).astype(int) < 2003
        ].index,
        inplace=True,
    )
    # Create escape character for sql request
    volcanoes.replace("'", "''", regex=True, inplace=True)
    # replace coma by points in number
    volcanoes.replace(",", ".", regex=True, inplace=True)

    volcanoes = pd.DataFrame(
        {
            "Volcano Number": volcanoes['Volcano Number'].to_numpy(),
            "Volcano Name": volcanoes['Volcano Name'].to_numpy(),
            "Country": volcanoes['Country'].to_numpy(),
            "Last Known Eruption": volcanoes['Last Known Eruption'].to_numpy(),
            "Latitude": volcanoes['Latitude'].to_numpy(),
            "Longitude": volcanoes['Longitude'].to_numpy(),
            "Elevation (m)": volcanoes['Elevation (m)'].to_numpy(),
        },
    )

    return volcanoes



base_path = '.'
folder_name = ''
run_date = ''
conf_file = ''
fig_name_detections = ''
fig_name_eruptions = ''
date_start = ''
dete_end = ''

date_start = datetime(2011, 1, 1).replace(tzinfo=timezone.utc)
date_end = datetime(2012, 1, 1).replace(tzinfo=timezone.utc)

run_date = sys.argv[1]
conf_file = [x for x in filter(lambda x: x.split('.')[-1] == 'toml',
                   os.listdir(join(base_path, 'results', sys.argv[1])))][0]

fig_name = f"{run_date}-results.png"

detections = pd.read_pickle(
    join(
        base_path,
        "compiled_data/detections.pkl"
        )
    )

cfg_stats = pd.read_csv('../cfg/stations.csv')

cfg_volcs = get_volcanoes_from_gvp_database('../cfg/volcanoes.csv')

cfg_file = toml.load(join(
    base_path,
    "results",
    run_date,
    conf_file
    ))

ips = pd.read_pickle(
    join(
        base_path,
        f"results/{run_date}/",
        'ip_results.pkl'
        )
    )
ips = ips[ips['Datetime (UTC)'] <= date_end]
ips = ips[ips['Datetime (UTC)'] >= date_start]

eruptions = pd.read_pickle(
    join(
        base_path,
        f"results/{run_date}/",
        'eruption_results.pkl'
        )
    )
eruptions = eruptions[eruptions['Start Date (UTC)'] >= date_start]
eruptions = eruptions[eruptions['End Date (UTC)'] <= date_end]

assoc_sta_er = pd.read_pickle(
    join(
        base_path,
        f"results/{run_date}/assoc_sta_er.pkl"
        )
    )

# --------------------------------------------------------------------------
# Get stations and volcano from configuration file
# --------------------------------------------------------------------------
stat_list = cfg_file['STATIONS']['StationList']

# Sort by distance ---------------------------------------------------------
dist_list = []
if len(cfg_file['VOLCANOES']['VolcanoesList']) > 0:
    volcano_name = cfg_file['VOLCANOES']['VolcanoesList'][0]  # NOTE: should be 1
    if type(volcano_name) is int:
        volcano_name = cfg_volcs[cfg_volcs['Volcano Number'] == volcano_name]['Volcano Name'].values[0]

    for stat_i in stat_list:
        sta_lat = cfg_stats[cfg_stats['Station Name'] == stat_i]['Latitude'].values[0]
        sta_lon = cfg_stats[cfg_stats['Station Name'] == stat_i]['Longitude'].values[0]
        vol_lat = cfg_volcs[cfg_volcs['Volcano Name'] == volcano_name]['Latitude'].values[0]
        vol_lon = cfg_volcs[cfg_volcs['Volcano Name'] == volcano_name]['Longitude'].values[0]
        d, az, baz = gps2dist_azimuth(vol_lat, vol_lon, sta_lat, sta_lon)
        dist_list.append(d/1000)
elif len(cfg_file['VOLCANOES']['Regions']) > 0:
    lats = []
    lons = []
    # NOTE: Just one region assumed
    for volcano_name in cfg_file['VOLCANOES']['Regions'][0]:
        lats.append(cfg_volcs[cfg_volcs['Volcano Name'] == volcano_name]['Latitude'].values[0])
        lons.append(cfg_volcs[cfg_volcs['Volcano Name'] == volcano_name]['Longitude'].values[0])

    vol_lat = np.mean(lats)
    vol_lon = np.mean(lons)
    for stat_i in stat_list:
        sta_lat = cfg_stats[cfg_stats['Station Name'] == stat_i]['Latitude'].values
        sta_lon = cfg_stats[cfg_stats['Station Name'] == stat_i]['Longitude'].values
        d, az, baz = gps2dist_azimuth(vol_lat, vol_lon, sta_lat, sta_lon)
        dist_list.append(d/1000)

dist_list = np.asarray(dist_list)
stat_list = np.asarray(stat_list)
sorted_by_dist_ind = list(reversed(np.argsort(dist_list)))
sorted_dist_list = dist_list[sorted_by_dist_ind]
sorted_stat_list = stat_list[sorted_by_dist_ind]

sorted_dist_list_aux = []
sorted_stat_list_aux = []
for ax_i, sta in enumerate(sorted_stat_list):
    ip = ips[ips['Station Name'] == sta]
    if len(ip) == 0:
        ip = ip[ip['Station Name'] == get_stat_name(sta)]
    ip = ip[ip['IP'] >= 1]
    if len(ip) > 0:
        sorted_dist_list_aux.append(sorted_dist_list[ax_i])
        sorted_stat_list_aux.append(sorted_stat_list[ax_i])
num_stations = len(sorted_stat_list_aux)
sorted_dist_list = sorted_dist_list_aux
sorted_stat_list = sorted_stat_list_aux


# --------------------------------------------------------------------------
#                                  Plot
# --------------------------------------------------------------------------

# Font sizes ---------------------------------------------------------------
plt.rc("font", size=10)  # controls default text sizes
plt.rc("axes", titlesize=10)  # fontsize of the axes title
plt.rc("axes", labelsize=10)  # fontsize of the x and y labels
plt.rc("xtick", labelsize=10)  # fontsize of the tick labels
plt.rc("ytick", labelsize=10)  # fontsize of the tick labels
plt.rc("legend", fontsize=10)  # legend fontsize
plt.rc("figure", titlesize=12)  # fontsize of the figure title

# Plot labels --------------------------------------------------------------
subplot_order = ['(a)', '(b)', '(c)', '(d)', '(e)', '(f)', '(g)', '(h)',
                     '(i)', '(j)', '(k)']  # up to 12 stations, modify if needed
subplot_order = subplot_order[0:num_stations+2]

# Fig size and shape -------------------------------------------------------
fig_size_x = 8
fig_size_y = 7

fig, axes = plt.subplots(
    num_stations+1,
    3,
    figsize=(fig_size_x, fig_size_y),
    gridspec_kw = {"width_ratios" : [0.01, 1, 0.01]},
    sharex=True
)

# Start plot ---------------------------------------------------------------
for ax_i, sta in enumerate(sorted_stat_list):
    ip = ips[ips['Station Name'] == sta]
    if len(ip) == 0:
        ip = ip[ip['Station Name'] == get_stat_name(sta)]
    ip = ip[ip['IP'] >= 1]

    ax_ip  = axes[ax_i,1]
    if ax_i<num_stations:
        ax_ip.set_xticklabels([])
    ax_ip.spines[['right', 'top', 'bottom']].set_visible(False)
    ax_ip.grid(axis='x')

    amp_1km = ip['Average Source Amplitude (Pa)'].values
    date_times_amp = ip['Datetime (UTC)'].values

    # Plot IP data
    scale_pix = .1  # NOTE: check this parameter
    h1 = ax_ip.scatter(
        ip['Datetime (UTC)'],
        np.log10(ip['IP']),
        s           = ip['Number of Detections']/scale_pix,
        c           = ip['Mean Frequency (Hz)'],
        cmap        = plt.get_cmap('viridis'),
        marker      = '.',
        edgecolor   = 'none',
        alpha       = 0.8,
        zorder      = 1000,
        vmin        = 1,
        vmax        = 3,
    )

    dist = sorted_dist_list[ax_i]

    # add text as title
    ax_ip.annotate(f"{subplot_order[ax_i]} {sta} ({dist:.1f} km)",
                     (0.005   , 1.03), xycoords='axes fraction',
                    ha='left', va='bottom', weight='bold',
                    bbox=dict(boxstyle='square', facecolor='white', alpha=0, pad=0.0,
                    edgecolor='white',
                    zorder=100)
    )

# Plot Eruption data
ax_er = axes[num_stations, 1]
ax_er.spines[['right', 'top', 'bottom']].set_visible(False)
ax_er.grid(axis='x')

num_stats_det = {}
color_levels = ['yellowgreen', 'coral', 'red', 'purple']
for i, er_i in eruptions.iterrows():
    start_date_i = er_i['Start Date (UTC)']
    end_date_i = er_i['End Date (UTC)']
    estimated_amp_i = er_i['Estimated Amplitude [Pa]']
    # Get how many stations detected this eruption
    eruption_code = er_i['Eruption Code']
    assoc_stats = assoc_sta_er[assoc_sta_er['Eruption Code'] == eruption_code]
    assoc_stats = assoc_stats[assoc_stats['Detecting'] == 1]

    ax_er.axvspan(start_date_i, end_date_i, color=color_levels[len(assoc_stats['Station Name'].values)-1])
    ax_er.plot([start_date_i, end_date_i],
                 [estimated_amp_i, estimated_amp_i],
                 linewidth=2, color='black')

    for sta_i in assoc_stats['Station Name'].values:
        if sta_i not in num_stats_det:
            num_stats_det[sta_i] = [[start_date_i, end_date_i]]
        else:
            num_stats_det[sta_i].append([start_date_i, end_date_i])

for ax_i, sta in enumerate(sorted_stat_list):
    erup_times = num_stats_det[sta]
    for start_i, end_i in erup_times:
        axes[ax_i, 1].axvspan(start_i, end_i, color='gray', alpha=0.5)

# add text as title
ax_er.annotate(f"{subplot_order[num_stations]} VIS eruptive periods and maximum amplitude at 1 km",
                 (0.005   , 1.03), xycoords='axes fraction',
                ha='left', va='bottom', weight='bold',
                bbox=dict(boxstyle='square', facecolor='white', alpha=0, pad=0.0,
                edgecolor='white',
                zorder=100)
)

# Axes on the left
gs = axes[0,0].get_gridspec()
for ax_i in axes[0:num_stations,0]:
    ax_i.remove()
ax_ylabel = fig.add_subplot(gs[0:num_stations,0])
ax_ylabel.axis('off')
ax_ylabel.annotate("log(IP)", (0.5, 0.5),
    xycoords='axes fraction',
    ha='center', va='center', rotation=90, fontsize=12)

axes[num_stations,0].remove()

ax_ylabel_2 = fig.add_subplot(gs[num_stations,0])
ax_ylabel_2.annotate("Amp. [Pa]", (0.5, 0.5),
                   xycoords='axes fraction',
                   ha='center', va='center', rotation=90, fontsize=12)
ax_ylabel_2.axis('off')

axes[num_stations,2].remove()

# Colorbar on the right
gs = axes[1,0].get_gridspec()
for ax_i in axes[0:num_stations,2]:
    ax_i.remove()
cax_pmcc = fig.add_subplot(gs[0:num_stations,2])
cbar = fig.colorbar(h1, cax=cax_pmcc)
cbar.ax.set_ylabel('Mean freq. [Hz]', fontsize=12)
cbar.ax.yaxis.set_major_formatter('{x:.1f}')

# Colorbar on the right for the eruptive periods
ax_er = fig.add_subplot(gs[num_stations,2])
ax_er.axis('off')
rect1 = plt.Rectangle((0.08, 0.05), 0.8, 0.3, facecolor=color_levels[0])
rect2 = plt.Rectangle((0.08, 0.35), 0.8, 0.3, facecolor=color_levels[1])
rect3 = plt.Rectangle((0.08, 0.65), 0.8, 0.3, facecolor=color_levels[2])
rect4 = plt.Rectangle((0.0, 0.05), 1., 0.9, facecolor=None, edgecolor='black', linewidth=2, fill=False, zorder=0)
ax_er.add_patch(rect1)
ax_er.add_patch(rect2)
ax_er.add_patch(rect3)
ax_er.add_patch(rect4)

ax_er.annotate('1', (0.6, 0.2), (2.5, 0.2), ha='left', va='center', arrowprops=dict(arrowstyle="-"))

ax_er.annotate('2', (0.6, 0.5), (2.5, 0.5), ha='left', va='center', arrowprops=dict(arrowstyle="-"))

ax_er.annotate('3', (0.6, 0.8), (2.5, 0.8), ha='left', va='center', arrowprops=dict(arrowstyle="-"))
ax_er.text(7, 0.5, '# stations', rotation=90, va='center', fontsize=12)

axes[num_stations,1].xaxis.set_major_formatter(
    mdates.ConciseDateFormatter(axes[num_stations,1].xaxis.get_major_locator()))
plt.subplots_adjust(hspace=0.7, wspace=0)

# Save the plot
plt.tight_layout()
fig_name = join('figures', fig_name)
plt.savefig(fig_name, dpi=300)
print(f"-> Figure <{fig_name}> saved.")
