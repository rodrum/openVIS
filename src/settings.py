"""
App settings
"""

import os
from os.path import dirname, abspath, join, exists
import sys
from datetime import datetime
import shutil
import toml

from src.logger import logger

# Build paths inside the project like this: os.path.join(BASE_DIR, ...)
BASE_DIR = dirname(dirname(abspath(__file__)))

# ==========================================================
# Configuration file
config = toml.load(join(BASE_DIR, 'cfg', 'vis_config.toml'))
# ==========================================================

# Start and end datetimes for analysis
START_DATE = config['DATES']['StartDate']
END_DATE = config['DATES']['EndDate']

# Volcanoes selected
VOLCANO_TABLE = join(BASE_DIR, 'cfg', config['VOLCANOES']['VolcanoesTable'])
MONITORING_AREA = config['VOLCANOES']['Area']
REGIONS = config['VOLCANOES']['Regions']
VOLCANOES = config['VOLCANOES']['VolcanoesList']

# Stations selected
STATIONS_TABLE = join(BASE_DIR, 'cfg', config['STATIONS']['StationsTable'])
STATIONS = [
    sta.upper() for sta in config['STATIONS']['StationList']
]

# Paths -----------------------------------------------------------------------
BULL_PATH = config['PATHS']['Bulletins']
VEFF_PATH = config['PATHS']['VeffRatios']
BAZDEV_PATH = config['PATHS']['BackAziInterp']
DATA_PATH = config['PATHS']['CompiledData']
RESULTS_PATH = config['PATHS']['Results']
list_of_paths = [BULL_PATH, VEFF_PATH, BAZDEV_PATH, DATA_PATH, RESULTS_PATH]
# Check they exist
for path_i in list_of_paths: 
    if path_i is False:
        continue
    elif not exists(path_i):
        logger.warning(f"Can't find {path_i}.")
        #sys.exit(0)
        os.makedirs(path_i, exist_ok=True)
        logger.warning(f"Created {path_i}.")
    else:
        logger.info(f"Path {path_i} found")


# 'Smart' results subfolder (If only)
datetime_now_str = datetime.now().strftime('%Y%m%dT%H%M%S')
this_results_folder = join(RESULTS_PATH, datetime_now_str)
if not exists(this_results_folder):
    os.mkdir(this_results_folder)
    logger.info(f"Folder <{this_results_folder}> will store your results"
                " and the configuration file")
    RESULTS_PATH = this_results_folder
else:
    logger.error(f"Path {this_results_folder} exists!")
    sys.exit(0)

# Formats ---------------------------------------------------------------------

VEFF_FORMAT = config['FORMATS']['VeffFormat']
allowed_veff_formats = ['BGR', 'CLIM']
if VEFF_FORMAT is False:
    logger.warning(f"VeffFormat if False -> veff-ratio = 1")
elif VEFF_FORMAT not in allowed_veff_formats:
    logger.error(f"VeffFormat not in: {allowed_veff_formats}")
    sys.exit(0)


# Processing parameters -------------------------------------------------------
# Minimum time interval between 2 notifications
NOTIFICATION_INTERVAL = config['PROCESSING']['NotificationInterval']
# Minimal IP value to trigger a notification
IP_THRESHOLD = config['PROCESSING']['IPThreshold']
# Maximal value of amplitude for a detection to be considered trustworthy
MAX_AMP = config['PROCESSING']['MaxAmp']
# Time interval between 2 IP calculations
IP_TIME_INTERVAL = config['PROCESSING']['WindowLength']
# Time window for a VIS analysis
ANALYSIS_TIME_INTERVAL = config['PROCESSING']['Timeshift']
# Delay without detection to consider an eruption as 'ended'
DELTA_CLOSE = config['PROCESSING']['DeltaClose']
# Maximal distance volcano to station
MAX_DIST = config['PROCESSING']['MaxDist']
# Azimuthal deviation tolerance
DAZIM = config['PROCESSING']['Dazim']
# Force DAZIM instead of daily dependent values
FORCE_DAZIM = config['PROCESSING']['ForceDazim']
# Number of STD widths if using back-azimuth estimation derived tolerance
NUM_BAZ_STD = None
if BAZDEV_PATH is not False:
    NUM_BAZ_STD = config['PROCESSING']['NumBazStd']
# Persistency threshold in percentage
PERSISTENCY_THRESHOLD = config['PROCESSING']['PersistencyThreshold']
# Max distance at which an eruption is generated with max confidence level
STA_CLOSE_CONFIDENCE = config['PROCESSING']['StaCloseConfidence']
# Minima signal central frequency
MIN_MEAN_FREQ = config['PROCESSING']['MinMeanFreq']
# Max signal mean frequency
MAX_MEAN_FREQ = config['PROCESSING']['MaxMeanFreq']
# Reference propagation speed of sound
REF_SPEED = config['PROCESSING']['RefSpeed']

# ============================================
# Save this configuration for future reference
# ============================================
num_volcs = 0
num_regions = 0
num_volcs = len(VOLCANOES)
if REGIONS != []:
    num_regions = len(REGIONS)

num_stats = len(STATIONS)

start_date_str = datetime.strftime(START_DATE, '%Y%m%dT%H%M%S')
end_date_str = datetime.strftime(END_DATE, '%Y%m%dT%H%M%S')

veff_ratio = 'Noveff'
if VEFF_PATH is not False:
    if VEFF_FORMAT == 'BGR':
        veff_ratio = 'BGRveff'
    elif VEFF_FORMAT == 'CLIM':
        veff_ratio = 'CLIMveff'

bazdev = 'Nobazdev'
if BAZDEV_PATH is not False:
    bazdev = 'UsingBazdev'

select_volcs = ''
if MONITORING_AREA != []:
    select_volcs = ".Area"  # NOTE: could use better descriptior
else:
    select_volcs = f".{num_volcs}volcs.{num_regions}regions"

name_config = f"{start_date_str}_{end_date_str}" +\
            f"{select_volcs}" +\
            f".{num_stats}stats.{veff_ratio}.{bazdev}.toml"

# Save file in DATA_PATH for future reference
shutil.copyfile(join(BASE_DIR, 'cfg', 'vis_config.toml'),
                join(BASE_DIR, this_results_folder, name_config))
