from unittest import TestCase
import tempfile
import os
import h5py

import numpy as np
from numpy.testing import assert_allclose, assert_equal, assert_array_equal

from km3net_testdata import data_path
import km3io as ki

import km3astro.io as io

# TestObject = ki.OfflineReader(data_path("offline/km3net_offline.root"))


class TestIO(TestCase):
    def setUp(self):
        self.root_file_path = data_path(
            "offline/mcv5.1.genhen_anumuNC.sirene.jte.jchain.aashower.sample.root"
        )
        self.TestObject = ki.OfflineReader(self.root_file_path)

        self.is_mc_data = True

        self.tmpdirname = tempfile.TemporaryDirectory()
        self.output_file_path = os.path.join(self.tmpdirname.name, "test_output.h5")
        self.hdf5_file_path = data_path(
            "hdf5/mcv5.1.genhen_anumuNC.sirene.jte.jchain.aashower.sample.hdf5"
        )

    def tearDown(self):
        # Cleanup: close the temporary directory
        self.tmpdirname.cleanup()

    def test_dir_to_spherical(self):
        # theta_exp = np.array([2.63188097, 2.55929208, 2.99414125, 2.88690687, 2.53314835, 2.21737661, 2.12586832, 2.17786349, 2.93141631, 2.4034084])
        theta_exp = np.array(
            [
                1.67854698,
                1.87939,
                1.6143014,
                0.73193345,
                1.59645405,
                0.5448071,
                1.60235583,
                0.82485171,
                1.16691096,
                1.74582966,
            ]
        )
        # phi_exp = np.array([4.788028156,  1.8949051577, 4.1228998863, 4.2939973467, 1.2371240233, 5.25252562,   2.141405076,  3.6462426279, 0.9341179294, 1.3128858853])
        phi_exp = np.array(
            [
                5.90517876,
                5.56053997,
                4.69128142,
                5.70289304,
                4.55218511,
                5.81203723,
                1.4096874,
                3.62824618,
                4.67115027,
                1.48297615,
            ]
        )
        theta, phi = io.dir_to_spherical(ki.tools.best_jmuon(self.TestObject.tracks))

        assert_allclose(theta, theta_exp)
        assert_allclose(phi, phi_exp)

    def test_get_utc_timeslice(self):
        # times_exp = np.array([1567036818.2 , 1567036818.3 , 1567036820.2 , 1567036816.5 , 1567036816.5 , 1567036816.5 , 1567036822.2 , 1567036818.5 , 1567036818.5 , 1567036820.4])
        times_exp = np.array(
            [175.5, 199.5, 295.6, 339.3, 560.1, 572.9, 574.6, 589.2, 627.6, 687.6]
        )
        times, astropytimes = io.get_utc_timeslice(self.TestObject)

        assert_allclose(times, times_exp)
        assert_allclose(np.array(astropytimes.value), times_exp)

    def test_get_utc_tracktime(self):
        utc_tracktimes_exp = np.array(
            [
                1567036818.27,
                1567036818.38,
                1567036820.22,
                1567036816.52,
                1567036816.55,
                1567036816.58,
                1567036822.25,
                1567036818.55,
                1567036818.6,
                1567036820.45,
            ]
        )
        timeslice_exp = np.array(
            [
                1567036818.2,
                1567036818.3,
                1567036820.2,
                1567036816.5,
                1567036816.5,
                1567036816.5,
                1567036822.2,
                1567036818.5,
                1567036818.5,
                1567036820.4,
            ]
        )
        tracktimes = ki.tools.best_jmuon(self.TestObject.tracks).t.to_numpy()

        utc_tracktimes, astropy_utc_tracktimes = io.get_utc_tracktime(
            timeslice_exp, tracktimes
        )

        assert_allclose(utc_tracktimes, utc_tracktimes_exp)
        assert_allclose(np.array(astropy_utc_tracktimes.value), utc_tracktimes_exp)

    def test_dump_header(self):
        header = io.dump_header(self.TestObject)
        self.assertIsNotNone(header)

    def test_generate_id_table(self):
        # timeslices_exp = np.array([1567036818.2 , 1567036818.3 , 1567036820.2 , 1567036816.5 , 1567036816.5 , 1567036816.5 , 1567036822.2 , 1567036818.5 , 1567036818.5 , 1567036820.4])
        timeslices_exp = np.array(
            [175.5, 199.5, 295.6, 339.3, 560.1, 572.9, 574.6, 589.2, 627.6, 687.6]
        )
        evt_ids_exp = np.array(range(0, 10))

        id_table, evt_ids, timeslices = io.generate_id_table(self.TestObject)

        dtype = id_table.dtype
        ncols = len(id_table.colnames)
        nrows = len(id_table)

        dtype_exp = np.dtype([("event_id", "<i8"), ("timeslice_utc_time", "O")])
        ncols_exp = 2
        nrows_exp = 10

        assert_array_equal(dtype, dtype_exp)

        assert_allclose(timeslices, timeslices_exp)
        assert_allclose(evt_ids, evt_ids_exp)

        assert_allclose(ncols, ncols_exp)
        assert_allclose(nrows, nrows_exp)

    def test_check_for_mc(self):
        is_mc_data = io.check_for_mc(self.TestObject)
        is_mc_data_exp = True
        assert_equal(is_mc_data, is_mc_data_exp)

    def test_generate_mc_table(self):
        evt_ids_exp = np.array(range(0, 10))

        if self.is_mc_data == True:
            mc_table = io.generate_mc_table(self.TestObject, evt_ids_exp)

            dtype = mc_table.dtype
            ncols = len(mc_table.colnames)
            nrows = len(mc_table)

            dtype_exp = np.dtype(
                [
                    ("event_id", "<i8"),
                    ("energy", "<f8"),
                    ("pdg_id", "<i4"),
                    ("normalized_weight", "<f8"),
                    ("theta_detectorframe", "<f8"),
                    ("phi_detectorframe", "<f8"),
                    ("pos_x", "<f8"),
                    ("pos_y", "<f8"),
                    ("pos_z", "<f8"),
                ]
            )
            ncols_exp = 9
            nrows_exp = 10

            assert_allclose(ncols, ncols_exp)
            assert_allclose(nrows, nrows_exp)
            assert_array_equal(dtype, dtype_exp)

            # Check shape and datatype of id_table
        else:
            print("self.is_mc_data == False")

    def test_generate_reco_table(self):
        timeslices_exp = np.array(
            [
                1567036818.2,
                1567036818.3,
                1567036820.2,
                1567036816.5,
                1567036816.5,
                1567036816.5,
                1567036822.2,
                1567036818.5,
                1567036818.5,
                1567036820.4,
            ]
        )
        evt_ids_exp = np.array(range(0, 10))

        reco_table = io.generate_reco_table(
            self.TestObject,
            evt_ids_exp,
            timeslices_exp,
            reco_types=["jmuon", "aashower", "jshower", "dusjshower"],
        )

        dtype = reco_table.dtype
        ncols = len(reco_table.colnames)
        nrows = len(reco_table)

        dtype_exp = np.dtype(
            [
                ("event_id", "<i8"),
                ("rec_type", "<i8"),
                ("rec_stage", "<i8"),
                ("tracktime_utc", "O"),
                ("tracktime_ns", "<f8"),
                ("energy", "<f8"),
                ("likelihood", "<f8"),
                ("theta_detectorframe", "<f8"),
                ("phi_detectorframe", "<f8"),
                ("pos_x", "<f8"),
                ("pos_y", "<f8"),
                ("pos_z", "<f8"),
            ]
        )
        ncols_exp = 12
        nrows_exp = 10

        assert_allclose(ncols, ncols_exp)
        assert_allclose(nrows, nrows_exp)
        assert_array_equal(dtype, dtype_exp)

    def test_root_to_hdf5(self):

        io.root_to_hdf5(self.root_file_path, output_file=self.output_file_path)

        # Check if the HDF5 file is created
        self.assertTrue(os.path.exists(self.output_file_path))

        with h5py.File(self.output_file_path, "r") as h5file:
            # Check if required groups and datasets exist
            self.assertTrue("RECO" in h5file)
            self.assertTrue("HEADER" in h5file)
            self.assertTrue("ID" in h5file)
            self.assertTrue("RECO_EVENTS" in h5file["RECO"])
            self.assertTrue("MC_EVENTS" in h5file["MC"])

    def test_load_hdf5_tables(self):

        tables = io.load_hdf5_tables(self.hdf5_file_path)

        # Check if the tables are loaded properly
        self.assertIsNotNone(tables.header_table)
        self.assertIsNotNone(tables.id_table)
        self.assertIsNotNone(tables.reco_table)
        # Check for optional table
        self.assertIsNotNone(tables.mc_table)

        # Expected shape information
        expected_id_table_num_rows = 10
        expected_id_table_colnames = ["event_id", "timeslice_utc_time"]

        expected_mc_table_num_rows = 10
        expected_mc_table_colnames = [
            "event_id",
            "energy",
            "pdg_id",
            "normalized_weight",
            "theta_detectorframe",
            "phi_detectorframe",
            "pos_x",
            "pos_y",
            "pos_z",
        ]

        expected_reco_table_num_rows = 10
        expected_reco_table_colnames = [
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
        ]

        # Check id_table
        self.assertEqual(
            len(tables.id_table),
            expected_id_table_num_rows,
            f"ID table row count mismatch. Expected {expected_id_table_num_rows}, got {len(tables.id_table)}",
        )
        self.assertListEqual(
            tables.id_table.colnames,
            expected_id_table_colnames,
            f"ID table column names mismatch. Expected {expected_id_table_colnames}, got {tables.id_table.colnames}",
        )

        # Check mc_table
        if tables.mc_table is not None:
            self.assertEqual(
                len(tables.mc_table),
                expected_mc_table_num_rows,
                f"MC table row count mismatch. Expected {expected_mc_table_num_rows}, got {len(tables.mc_table)}",
            )
            self.assertListEqual(
                tables.mc_table.colnames,
                expected_mc_table_colnames,
                f"MC table column names mismatch. Expected {expected_mc_table_colnames}, got {tables.mc_table.colnames}",
            )

        # Check reco_table
        self.assertEqual(
            len(tables.reco_table),
            expected_reco_table_num_rows,
            f"RECO table row count mismatch. Expected {expected_reco_table_num_rows}, got {len(tables.reco_table)}",
        )
        self.assertListEqual(
            tables.reco_table.colnames,
            expected_reco_table_colnames,
            f"RECO table column names mismatch. Expected {expected_reco_table_colnames}, got {tables.reco_table.colnames}",
        )
