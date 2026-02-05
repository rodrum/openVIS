"""
This module details the models used for the database.
"""

import sys
import pandas as pd
from datetime import datetime, timezone
from os.path import join

sys.path.append(".")

from src.settings import RESULTS_PATH
from src.logger import logger

def create_output_dataframes():
    """FILLME"""
    eruption_results = pd.DataFrame(
        {
            'Eruption Code': ['None'],
            'Start Date (UTC)': [datetime(1, 1, 1, tzinfo=timezone.utc)],
            'End Date (UTC)': [datetime(1, 1, 1, tzinfo=timezone.utc)],
            'Detection Date (UTC)': [datetime(1, 1, 1, tzinfo=timezone.utc)],
            'Confidence Level': [0],
            'Last Notification (UTC)': [datetime(1, 1, 1, tzinfo=timezone.utc)],
            'Revision': [0],
            'Status': ['None'],
            'Estimated Amplitude [Pa]': [0.0],
        }
    )

    assoc_sta_er = pd.DataFrame(
        {
            'Station Name': ['None'],
            'Eruption Code': ['None'],
            'Num. Detections': [0],
            'Max. Amp. [Pa]': [0.0],
            'Estimated Amp. [Pa]': [0.0],
            'Detecting': [0]
        }
    )

    assoc_volc_er = pd.DataFrame(
        {
            'Eruption Code': ['None'],
            'Volcano Code': [0],
            'Volcano Name': ['None']
        }
    )

    ip_results = pd.DataFrame(
        {
            'Station Name': ['None'],
            'Volcano Code': [0],
            'IP': [0.0],
            'Mean Amplitude (Pa)': [0.0],
            'Number of Detections': [0],
            'Veff-ratio': [0.0],
            'Back-Azimuth (deg)': [0.0],
            'Distance (km)': [0.0],
            'Datetime (UTC)': [datetime(1, 1, 1, tzinfo=timezone.utc)],
            'Mean Source Amplitude (Pa)': [0.0],
            'Persistency': [0.0],
            'Mean Frequency (Hz)': [0.0]
        }
    )

    eruption_results_file = join(RESULTS_PATH, 'eruption_results.pkl')
    ip_results_file = join(RESULTS_PATH, 'ip_results.pkl')
    assoc_sta_er_file = join(RESULTS_PATH, 'assoc_sta_er.pkl')
    assoc_volc_er_file = join(RESULTS_PATH, 'assoc_volc_er.pkl')

    eruption_results.to_pickle(eruption_results_file)
    logger.info(f"Created {eruption_results_file}")
    ip_results.to_pickle(ip_results_file)
    logger.info(f"Created {ip_results_file}")
    assoc_sta_er.to_pickle(assoc_sta_er_file)
    logger.info(f"Created {assoc_sta_er_file}")
    assoc_volc_er.to_pickle(assoc_volc_er_file)
    logger.info(f"Created {assoc_volc_er_file}")


if __name__ == '__main__':
    create_output_dataframes()
