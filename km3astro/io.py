"""
io.py
Converting km3net root files to astropy tables in hdf5 files
High-level data relevant for astro analyses is extracted and saved

author: 
Hannes Warnhofer 
hannes.warnhofer@fau.de

"""
import os

import time
from h5py import File, Group

# import h5py

import awkward as ak
import numpy as np

# from numpy import array, full_like, transpose, uint64

from km3io import offline
from km3io import tools
import km3io.definitions as kd
import datetime
import astropy.units as u
from astropy.table import Table, vstack
from astropy.time import Time
from astropy.io.misc.hdf5 import write_table_hdf5, read_table_hdf5


def dir_to_spherical(tracks):
    """
    Convert directional coordinates (x, y, z) to spherical coordinates (theta, phi).

    Parameters
    ----------
    tracks : numpy.array
        Array containing track information, OfflineReader.tracks object

    Returns
    -------
    tuple
        Tuple containing spherical coordinates (theta, phi).
    """
    vec = np.array([tracks.dir_x, tracks.dir_y, tracks.dir_z])
    vec = np.transpose(vec)
    theta_val = tools.theta(vec)
    phi_val = tools.phi(vec)
    return theta_val, phi_val


def get_utc_timeslice(OfflineObject):

    """
    Calculate UTC time of the timeslices where the events were put based on input data.

    Parameters:
    ----------
    OfflineObject : km3io.offline.OfflineReader
        Input data. Offline data runs: km3io.offline.OfflineReader Object.

    Returns:
    ----------
    tuple
        Tuple containing:
        - times : list
            List of UTC timeslices in seconds
        - astropy_times : list
            List of astropy.time.Time objects of times
    """
    tsec = OfflineObject.t_sec
    tns = OfflineObject.t_ns
    S_TO_NS = np.uint64(1e9)

    # print("calculating timesclices ...")

    times = (
        tsec + tns / S_TO_NS
    )  # Combine seconds and nanoseconds to get time in seconds
    astropy_times = Time(tsec, tns / S_TO_NS, format="unix")

    return times, astropy_times


def get_utc_tracktime(timeslice, tracktimes):
    """
    Calculate UTC time of the track times based on timeslices and track times.

    Parameters:
    ----------
    timeslice : list
        List of timeslices in seconds.
    tracktimes : list
        List of track times in nanoseconds.

    Returns:
    ----------
    tuple
        Tuple containing:
        - times : list
            List of UTC track times in seconds.
        - astropy_times : list
            List of astropy.time.Time objects of track times.
    """
    S_TO_NS = np.uint64(1e9)

    # print("calculating tracktimes ...")

    times = (
        timeslice + tracktimes / S_TO_NS
    )  # Combine seconds and nanoseconds to get time in seconds
    astropy_times = Time(timeslice, tracktimes / S_TO_NS, format="unix")

    return times, astropy_times


def dump_header(OfflineObject):
    """
    Create an astropy table containing header information from a file.

    Parameters
    ----------
    OfflineObject : km3io.offline.OfflineReader
        OfflineReader object representing the input file.

    Returns
    -------
    astropy.table.Table
        Table containing header information converted to strings.
    """
    # Create an empty table
    header_table = Table(names=["Parameter", "Value"], dtype=[str, str])

    # Loop through each attribute of the header object and add it to the table
    for key in dir(OfflineObject.header):
        header_table.add_row([key, str(getattr(OfflineObject.header, key))])
    # print("header_table created ...")

    return header_table


def generate_id_table(OfflineObject, evt_ids_start=0):
    """
    Generate event IDs and Monte Carlo data astropy tables from the given file.

    Parameters
    ----------
    OfflineObject : km3io.offline.OfflineReader
        OfflineReader object representing the input file.
    evt_ids_start : int, optional
        Starting event ID.

    Returns
    -------
    tuple
        Tuple containing ID table, event IDs as numpy array, and timeslices.
        - ID table : astropy.table.Table
            Table containing Event_ID and Timeslice_UTC_Time columns.
        - event_ids_np : numpy.ndarray
            Array of event IDs.
        - timeslices : list
            List of timeslices in seconds.
    """
    # Setting up event ids
    event_ids = range(evt_ids_start, len(OfflineObject) + evt_ids_start)
    event_ids_np = np.array(event_ids)
    # print(f"{event_ids = }")

    # Getting utc time of the timeslice the events were triggered
    timeslices, astropy_timeslices = get_utc_timeslice(OfflineObject)

    id_table = Table(
        [event_ids_np, astropy_timeslices], names=["event_id", "timeslice_utc_time"]
    )
    # print("id_table created ...")

    return id_table, event_ids_np, timeslices


def generate_mc_table(OfflineObject, evt_ids):
    """
    Generate monte carlo data astropy table

    Parameters
    ----------
    OfflineObject: km3io.OfflineReader
        OfflineReader Object representing the input file.
    evt_ids: numpy.array
        Array containing event IDs.

    Returns
    -------
    astropy.table.Table
        astropy table containing mc data.
    """

    # Getting mc neutrino
    nu = OfflineObject.mc_trks[:, 0]
    # Getting energy and origin of the neutrino
    E = nu.E.to_numpy()
    pdgid = nu.pdgid.to_numpy()
    theta, phi = dir_to_spherical(nu)
    posx = nu.pos_x.to_numpy()
    posy = nu.pos_y.to_numpy()
    posz = nu.pos_z.to_numpy()
    weight = OfflineObject.w.to_numpy()[
        :, 1
    ]  # check km3net dataformat definitions w2list
    num_gen_events = OfflineObject.header.genvol.numberOfEvents
    norm_weight = weight / num_gen_events

    mc_table = Table(
        [
            evt_ids,
            E * u.GeV,
            pdgid,
            norm_weight,
            theta * u.rad,
            phi * u.rad,
            posx,
            posy,
            posz,
        ],
        names=[
            "event_id",
            "energy",
            "pdg_id",
            "normalized_weight",
            "theta_detectorframe",
            "phi_detectorframe",
            "pos_x",
            "pos_y",
            "pos_z",
        ],
    )
    # print("mc_table created ...")
    return mc_table


def generate_reco_table(OfflineObject, evt_ids, timeslices, reco_types):
    """
    Generate reconstruction data astropy tables for specified reconstruction types.

    Parameters
    ----------
    OfflineObject : km3io.offline.OfflineReader
        OfflineReader object representing the input file.
    evt_ids : numpy.array
        Array containing event IDs.
    timeslices : list
        List of timeslices in seconds.
    reco_types : list of str
        List of reconstruction types. Allowed types: 'jmuon','jshower','aashower','dusjshower'.

    Returns
    -------
    astropy.table.Table
        Combined reconstruction table.
    """
    reco_tables = []
    for reco_type in reco_types:
        # Fetching best fits for the current reconstruction type
        # print(f"dealing with {reco_type} reco_type ...")
        if reco_type == "jmuon":
            if (
                any(
                    element is not False
                    for element in tools.has_jmuon(OfflineObject.tracks)
                )
                == False
            ):
                print(f"No {reco_type} reconstructions for this file!")
                continue
            else:
                best_fit = tools.best_jmuon(OfflineObject.tracks)
                reco_type_id = (
                    kd.reconstruction.JPP_RECONSTRUCTION_TYPE
                )  #'JPP_RECONSTRUCTION_TYPE'
                reco_stage_id = kd.reconstruction.JMUONBEGIN  #'JMUONBEGIN'
                mask_reco = tools.has_jmuon(OfflineObject.tracks)

        elif reco_type == "aashower":
            if (
                any(
                    element is not False
                    for element in tools.has_aashower(OfflineObject.tracks)
                )
                == False
            ):
                print(f"No {reco_type} reconstructions for this file!")
                continue
            else:
                best_fit = tools.best_aashower(OfflineObject.tracks)
                reco_type_id = (
                    kd.reconstruction.AANET_RECONSTRUCTION_TYPE
                )  #'AANET_RECONSTRUCTION_TYPE'
                reco_stage_id = kd.reconstruction.AASHOWERBEGIN  #'AASHOWERBEGIN'
                mask_reco = tools.has_aashower(OfflineObject.tracks)

        elif reco_type == "jshower":
            if (
                any(
                    element is not False
                    for element in tools.has_jshower(OfflineObject.tracks)
                )
                == False
            ):
                print(f"No {reco_type} reconstructions for this file!")
                continue
            else:
                best_fit = tools.best_jshower(OfflineObject.tracks)
                reco_type_id = (
                    kd.reconstruction.JPP_RECONSTRUCTION_TYPE
                )  #'JPP_RECONSTRUCTION_TYPE'
                reco_stage_id = kd.reconstruction.JSHOWERBEGIN  #'JSHOWERBEGIN'
                mask_reco = tools.has_jshower(OfflineObject.tracks)

        elif reco_type == "dusjshower":
            if (
                any(
                    element is not False
                    for element in tools.has_dusjshower(OfflineObject.tracks)
                )
                == False
            ):
                print(f"No {reco_type} reconstructions for this file!")
                continue
            else:
                best_fit = tools.best_dusjshower(OfflineObject.tracks)
                reco_type_id = (
                    kd.reconstruction.DUSJ_RECONSTRUCTION_TYPE
                )  #'DUSJ_RECONSTRUCTION_TYPE'
                reco_stage_id = kd.reconstruction.DUSJSHOWERBEGIN  #'DUSJSHOWERBEGIN'
                mask_reco = tools.has_dusjshower(OfflineObject.tracks)

        # print(reco_type_id)

        best_fit = best_fit[mask_reco]
        event_ids_reco = evt_ids[mask_reco]

        # Extracting data from the best fits
        E = best_fit.E.to_numpy()
        lik = best_fit.lik.to_numpy()
        timedata = best_fit.t.to_numpy()  # ns
        theta, phi = dir_to_spherical(best_fit)
        posx = best_fit.pos_x.to_numpy()
        posy = best_fit.pos_y.to_numpy()
        posz = best_fit.pos_z.to_numpy()
        reco_type_id_reco = np.full_like(event_ids_reco, reco_type_id)
        reco_stage_id_reco = np.full_like(event_ids_reco, reco_stage_id)

        # Calculating UTC tracktimes
        utc_tracktimes, astropy_utc_tracktimes = get_utc_tracktime(
            np.array(timeslices)[mask_reco], timedata
        )

        # Creating table for current reconstruction type
        reco_table = Table(
            [
                event_ids_reco,
                reco_type_id_reco,
                reco_stage_id_reco,
                astropy_utc_tracktimes,
                timedata * u.ns,
                E * u.GeV,
                lik,
                theta * u.rad,
                phi * u.rad,
                posx,
                posy,
                posz,
            ],
            names=[
                "event_id",
                "rec_type",
                "rec_stage",
                "tracktime_utc",
                "tracktime_ns",
                "energy",
                "likelihood",
                "theta_detectorframe",
                "phi_detectorframe",
                "pos_x",
                "pos_y",
                "pos_z",
            ],
        )
        reco_tables.append(reco_table)

    reco_tables_combined = vstack(reco_tables)
    reco_tables_combined.sort("event_id")
    # print("reco_table created ...")
    return reco_tables_combined


def check_for_mc(OfflineObject):
    """
    Check if the file contains Monte Carlo (MC)-generated or real data based on the entries in mc_trks.id.

    Parameters
    ----------
    OfflineObject : km3io.offline.OfflineReader
        OfflineReader object representing the input file.

    Returns
    -------
    bool
        True if the file contains MC-generated data, False otherwise.
    """
    ids = OfflineObject.mc_trks.id
    if ak.any(ids):
        print("The file contains MC-generated data.")
        return True
    else:
        print(
            "The file contains real data. Or at least doesn't have any entries in mc_trks.id ..."
        )
        return False


def root_to_hdf5(
    file_path,
    reco_types=["jmuon", "aashower", "jshower", "dusjshower"],
    output_file=None,
    event_ids_start=0,
):
    """
    Load a KM3NeT offline file, process it to extract high-level data that is relevant for astrophysical analysis in various astropy tables, and store them in an HDF5 file.

    Parameters
    ----------
    file_path : str
        Path to the input file.
    reco_types : list of str, optional
        List of reconstruction types. Allowed types: 'jmuon','jshower','aashower','dusjshower'.
        Default: ['jmuon','aashower','jshower','dusjshower']
    h5path : str, optional
        Path for storing the HDF5 file. Default: "OutputTables/"
    event_ids_start : int, optional
        Starting event ID. Default: 0

    Returns
    -------
    None
    """
    start_time = time.time()

    # Accessing the filename of the root file in order to match the name of the hdf5 file that gets created
    folder_path, file_name = os.path.split(file_path)
    file_name = os.path.splitext(file_name)[0]
    if output_file == None:
        combined_table_path = file_name + ".h5"
    else:
        combined_table_path = output_file
    # Loading offline data with the km3io OfflineReader
    OfflineObject = offline.OfflineReader(file_path)

    # Checking if the data contains mc information that can be written or not
    is_mc_data = check_for_mc(OfflineObject)

    # Generating the individual tables
    header_table = dump_header(OfflineObject)
    print("header_table created...")
    id_table, event_ids, timeslices = generate_id_table(OfflineObject, event_ids_start)
    print("id_table created...")
    reco_table = generate_reco_table(OfflineObject, event_ids, timeslices, reco_types)
    print("reco_table created...")
    if is_mc_data == True:
        print("is_mc_data == True")
        mc_table = generate_mc_table(OfflineObject, event_ids)

    # Writing the tables to a hdf5 file
    with File(combined_table_path, "w") as h5file:
        reco_grp = h5file.create_group("RECO")
        write_table_hdf5(header_table, h5file, path="HEADER", serialize_meta=True)
        write_table_hdf5(id_table, h5file, path="ID", serialize_meta=True)
        write_table_hdf5(reco_table, reco_grp, path="RECO_EVENTS", serialize_meta=True)
        if is_mc_data == True:
            mc_grp = h5file.create_group("MC")
            write_table_hdf5(mc_table, mc_grp, path="MC_EVENTS", serialize_meta=True)
    print(f"File written: {combined_table_path}")

    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Total time: {elapsed_time:.3f} seconds")


def load_hdf5_tables(file_path):
    """
    Load an HDF5 file and assign tables to an object.
    The function expects a file that is generated with astropy_hdf5_from_root() and has the paths 'HEADER','ID','RECO_EVENTS' and optionally 'MC_EVENTS'

    Parameters
    ----------
    file_path : str
        Path to the HDF5 file.

    Returns
    -------
    object
        An object containing attributes header_table, id_table, reco_table, and optionally mc_table.
    """
    tables = type("Tables", (), {})()  # Creating an empty object

    try:
        with File(file_path, "r") as h5file:
            tables.header_table = read_table_hdf5(h5file, path="HEADER")
            tables.id_table = read_table_hdf5(h5file, path="ID")
            if "RECO/RECO_EVENTS" in h5file:
                print("RECO_EVENTS is present")
                tables.reco_table = read_table_hdf5(h5file, path="RECO/RECO_EVENTS")

            else:
                print("RECO_EVENTS is NOT present")
                # tables.reco_table = read_table_hdf5(h5file, path='RECO/RECO_EVENTS')
            if "MC/MC_EVENTS" in h5file:
                tables.mc_table = read_table_hdf5(h5file, path="MC/MC_EVENTS")
                print("MC_EVENTS is present")
            else:
                print("MC_EVENTS is NOT present")
                # tables.mc_table = read_table_hdf5(h5file, path='MC/MC_EVENTS')
                tables.mc_table = None
    except Exception as e:
        print(f"Error loading tables from {file_path}: {e}")

    return tables
