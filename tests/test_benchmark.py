import numpy as np
import pandas as pd

from km3net_testdata import data_path
from km3astro.coord import *

def test_coordinate_transformator(file0 = "", detector_ = "antares", detector_to_ = "antares", type_name = "coordinate_systems"):

    print("Starting Test Coordinate Transformator...")
    #take a benchmark and create new benchmark table

    if file0 == "":
        file0 = data_path("astro/antares_coordinate_systems_benchmark.csv")

    table_read = reader_from_file(file0)

    print("Read " + file0 + " file. \n")
    print(table_read)
    path = "/home/htedjditi/work/test_astro/"
    g_comment = "# benchmark table from https://git.km3net.de/km3py/km3net-testdata/-/tree/master/km3net_testdata/data/astro\n# values were obtained with the corresponding benchmark "

    has_object = False
    date_ = table_read['date'] 
    time_ = table_read['time']
    if 'object' in table_read.columns:
        object_ = table_read["object"]
        has_object = True

    if set(['phi','theta']).issubset(table_read.columns):

        table_loc_to_utm = transform_to_new_frame(table_read, "Detector", "UTM", detector_, detector_to_)
        table_loc_to_eq = transform_to_new_frame(table_read, "Detector", "equatorial", detector_, detector_to_)
        table_loc_to_gal = transform_to_new_frame(table_read, "Detector", "galactic", detector_, detector_to_)

        phi_ = table_read['phi']
        theta_ = table_read['theta']

        az_, ze_ = zip(*table_loc_to_utm.apply( lambda x: get_az_zenith(x.SkyCoord_new, detector_to_) , axis = 1))
        
        ra_, dec_ = zip(*table_loc_to_eq.apply( lambda x: get_ra_dec(x.SkyCoord_new) , axis = 1))     
        l_, b_ = zip(*table_loc_to_gal.apply( lambda x: get_l_b(x.SkyCoord_new) , axis = 1))     

        data_ = {'date': date_, 'time' : time_, 'phi' : phi_ , 'theta' : theta_ , 'azimuth' : az_ , 'zenith' : ze_ ,'RA-ICRS' : ra_ ,'DEC-ICRS' : dec_ ,'gal_lon' : l_ ,'gal_lat' : b_ , }
        if has_object == True:
            data_ = {'date': date_, 'time' : time_, 'object' : object_, 'phi' : phi_ , 'theta' : theta_ , 'azimuth' : az_ , 'zenith' : ze_ ,'RA-ICRS' : ra_ ,'DEC-ICRS' : dec_ ,'gal_lon' : l_ ,'gal_lat' : b_ , }
        
        df_ = pd.DataFrame(data=data_)
        
        comment = g_comment + "for 'phi' and 'theta'\n# other value are obtained with km3py_astro.\n# units: [deg]\n# theta, phi: track in horiz. UTM\n# zenith, azimuth: Source in horiz. UTM\n# DEC/RA: source in equatorial coos.\n# gal_lon, gal_lat: source in galactic coos.\n"
        name = path + detector_ + "_" + type_name + "_benchmark_km3py_from_track.csv"
        f = open(name, 'w')
        f.write(comment)
        
        df_.to_csv(f , index = False)
        f.close()

    if set(['azimuth','zenith']).issubset(table_read.columns):

        table_utm_to_loc = transform_to_new_frame(table_read, "UTM", "Detector", detector_, detector_to_)
        table_utm_to_eq = transform_to_new_frame(table_read, "UTM", "equatorial", detector_, detector_to_)
        table_utm_to_gal = transform_to_new_frame(table_read, "UTM", "galactic", detector_, detector_to_)

        az_ = table_read['azimuth']
        ze_ = table_read['zenith']

        phi_, theta_ = zip(*table_utm_to_loc.apply( lambda x: get_phi_theta(x.SkyCoord_new, detector_to_) , axis = 1))
        
        ra_, dec_ = zip(*table_utm_to_eq.apply( lambda x: get_ra_dec(x.SkyCoord_new) , axis = 1))     
        l_, b_ = zip(*table_utm_to_gal.apply( lambda x: get_l_b(x.SkyCoord_new) , axis = 1))     

        data_ = {'date': date_, 'time' : time_, 'phi' : phi_ , 'theta' : theta_ , 'azimuth' : az_ , 'zenith' : ze_ ,'RA-ICRS' : ra_ ,'DEC-ICRS' : dec_ ,'gal_lon' : l_ ,'gal_lat' : b_ , }
        if has_object == True:
            data_ = {'date': date_, 'time' : time_, 'object' : object_, 'phi' : phi_ , 'theta' : theta_ , 'azimuth' : az_ , 'zenith' : ze_ ,'RA-ICRS' : ra_ ,'DEC-ICRS' : dec_ ,'gal_lon' : l_ ,'gal_lat' : b_ , }
        
        df_ = pd.DataFrame(data=data_)
        
        comment = g_comment + "for 'azimuth' and 'zenith'\n# other value are obtained with km3py_astro.\n# units: [deg]\n# theta, phi: track in horiz. UTM\n# zenith, azimuth: Source in horiz. UTM\n# DEC/RA: source in equatorial coos.\n# gal_lon, gal_lat: source in galactic coos.\n"
        name = path + detector_ + "_" + type_name + "_benchmark_km3py_from_utm.csv"
        f = open(name, 'w')
        f.write(comment)
        
        df_.to_csv(f , index = False)
        f.close()

    if set(['RA-J2000','DEC-J2000']).issubset(table_read.columns):

        table_eq_to_utm = transform_to_new_frame(table_read, "equatorial", "UTM", detector_, detector_to_)
        table_eq_to_loc = transform_to_new_frame(table_read, "equatorial", "Detector", detector_, detector_to_)
        table_eq_to_gal = transform_to_new_frame(table_read, "equatorial", "galactic", detector_, detector_to_)

        ra_ = table_read['RA-J2000']
        dec_ = table_read['DEC-J2000']

        phi_, theta_ = zip(*table_eq_to_loc.apply( lambda x: get_phi_theta(x.SkyCoord_new, detector_to_) , axis = 1))

        az_, ze_ = zip(*table_eq_to_utm.apply( lambda x: get_az_zenith(x.SkyCoord_new, detector_to_) , axis = 1))
        
        l_, b_ = zip(*table_eq_to_gal.apply( lambda x: get_l_b(x.SkyCoord_new) , axis = 1))     

        data_ = {'date': date_, 'time' : time_, 'phi' : phi_ , 'theta' : theta_ , 'azimuth' : az_ , 'zenith' : ze_ ,'RA-J2000' : ra_ ,'DEC-J2000' : dec_ ,'gal_lon' : l_ ,'gal_lat' : b_ , }
        if has_object == True:
            data_ = {'date': date_, 'time' : time_, 'object' : object_, 'phi' : phi_ , 'theta' : theta_ , 'azimuth' : az_ , 'zenith' : ze_ ,'RA-J2000' : ra_ ,'DEC-J2000' : dec_ ,'gal_lon' : l_ ,'gal_lat' : b_ , }
        
        df_ = pd.DataFrame(data=data_)
        
        comment = g_comment + "for 'RA-J2000' and 'RA-J2000'\n# other value are obtained with km3py_astro.\n# units: [deg]\n# theta, phi: track in horiz. UTM\n# zenith, azimuth: Source in horiz. UTM\n# DEC/RA: source in equatorial coos.\n# gal_lon, gal_lat: source in galactic coos.\n"
        name = path + detector_ + "_" + type_name + "_benchmark_km3py_from_equatorial.csv"
        f = open(name, 'w')
        f.write(comment)
        
        df_.to_csv(f , index = False)
        f.close()

    if set(['gal_lon','gal_lat']).issubset(table_read.columns):

        table_gal_to_loc = transform_to_new_frame(table_read, "galactic", "Detector", detector_, detector_to_)
        table_gal_to_eq = transform_to_new_frame(table_read, "galactic", "equatorial", detector_, detector_to_)
        table_gal_to_utm = transform_to_new_frame(table_read, "galactic", "UTM", detector_, detector_to_)

        l_ = table_read['gal_lon']
        b_ = table_read['gal_lat']

        phi_, theta_ = zip(*table_utm_to_loc.apply( lambda x: get_phi_theta(x.SkyCoord_new, detector_to_) , axis = 1))
        az_, ze_ = zip(*table_eq_to_utm.apply( lambda x: get_az_zenith(x.SkyCoord_new, detector_to_) , axis = 1))

        ra_, dec_ = zip(*table_utm_to_eq.apply( lambda x: get_ra_dec(x.SkyCoord_new) , axis = 1))

        data_ = {'date': date_, 'time' : time_, 'phi' : phi_ , 'theta' : theta_ , 'azimuth' : az_ , 'zenith' : ze_ ,'RA-ICRS' : ra_ ,'DEC-ICRS' : dec_ ,'gal_lon' : l_ ,'gal_lat' : b_ , }
        if has_object == True:
            data_ = {'date': date_, 'time' : time_, 'object' : object_, 'phi' : phi_ , 'theta' : theta_ , 'azimuth' : az_ , 'zenith' : ze_ ,'RA-ICRS' : ra_ ,'DEC-ICRS' : dec_ ,'gal_lon' : l_ ,'gal_lat' : b_ , }
        
        df_ = pd.DataFrame(data=data_)
        
        comment = g_comment + "for 'gal_lon' and 'gal_lat'\n# other value are obtained with km3py_astro.\n# units: [deg]\n# theta, phi: track in horiz. UTM\n# zenith, azimuth: Source in horiz. UTM\n# DEC/RA: source in equatorial coos.\n# gal_lon, gal_lat: source in galactic coos.\n"
        name = path + detector_ + "_" + type_name + "_benchmark_km3py_from_galactic.csv"
        f = open(name, 'w')
        f.write(comment)
        
        df_.to_csv(f , index = False)
        f.close()

    print("... End of Test Coordinate Transformator")

    return 0

def test_Skycoord_separation(SC_true, SC_check):
    
    sep = SC_true.separation(SC_check)
    sep_deg = sep.deg
    return sep_deg

def test_benchmark_conversion(table_true, table_check):

    data_ = [table_true["SkyCoord_base"], table_check["SkyCoord_new"]]
    table_ = pd.concat(data_, axis = 1) 

    sep_table = table_.apply( lambda x: test_Skycoord_separation(x.SkyCoord_base, x.SkyCoord_new), axis = 1, result_type='expand')

    return sep_table

def test_angle_separation(file0 = "", detector_ = "antares", detector_to_ = "antares", to_print = False):

    print("Starting Angle Separation Test...")
    if file0 == "":
        file0 = data_path("astro/antares_coordinate_systems_benchmark.csv")

    print("Read " + file0 + " file. \n")

    table_read = reader_from_file(file0)
    print(table_read)

    angle_treshold = 0.02
    #Calculate value for all frame from initial benchmark value

    if set(['phi','theta']).issubset(table_read.columns):

        table_loc_to_utm = transform_to_new_frame(table_read, "Detector", "UTM", detector_, detector_to_)
        table_loc_to_eq = transform_to_new_frame(table_read, "Detector", "equatorial", detector_, detector_to_)
        table_loc_to_gal = transform_to_new_frame(table_read, "Detector", "galactic", detector_, detector_to_)

    if set(['azimuth','zenith']).issubset(table_read.columns):

        table_utm_to_loc = transform_to_new_frame(table_read, "UTM", "Detector", detector_, detector_to_)
        table_utm_to_eq = transform_to_new_frame(table_read, "UTM", "equatorial", detector_, detector_to_)
        table_utm_to_gal = transform_to_new_frame(table_read, "UTM", "galactic", detector_, detector_to_)

    if set(['RA-J2000','DEC-J2000']).issubset(table_read.columns):

        table_eq_to_utm = transform_to_new_frame(table_read, "equatorial", "UTM", detector_, detector_to_)
        table_eq_to_loc = transform_to_new_frame(table_read, "equatorial", "Detector", detector_, detector_to_)
        table_eq_to_gal = transform_to_new_frame(table_read, "equatorial", "galactic", detector_, detector_to_)

    if set(['gal_lon','gal_lat']).issubset(table_read.columns):

        table_gal_to_loc = transform_to_new_frame(table_read, "galactic", "Detector", detector_, detector_to_)
        table_gal_to_eq = transform_to_new_frame(table_read, "galactic", "equatorial", detector_, detector_to_)
        table_gal_to_utm = transform_to_new_frame(table_read, "galactic", "UTM", detector_, detector_to_)

    #Start the Angle separation calculation
    print("Checking Phi/Theta")
    
    if set(['phi','theta']).issubset(table_read.columns):
        
        if set(['azimuth','zenith']).issubset(table_read.columns):
            sep_utm_to_loc = test_benchmark_conversion(table_loc_to_utm, table_utm_to_loc)
            
            mean_ = sep_utm_to_loc.mean()
            min_ = sep_utm_to_loc.min()
            max_ = sep_utm_to_loc.max()
            
            if to_print == True:
                print("UTM TO LOC")
                print(sep_utm_to_loc)
                print("Mean = " + str(mean_) + " Max = " + str(max_) + " Min = " + str(min_))
            
            if max_ > angle_treshold:
                raise Exception("Error: Maximum angle separation = " + str(max_) + " > " + str(angle_treshold))


        if set(['RA-J2000','DEC-J2000']).issubset(table_read.columns):
            sep_eq_to_loc = test_benchmark_conversion(table_loc_to_eq, table_eq_to_loc)

            mean_ = sep_eq_to_loc.mean()
            min_ = sep_eq_to_loc.min()
            max_ = sep_eq_to_loc.max()

            if to_print == True:
                print("EQ TO LOC")
                print("Mean = " + str(mean_) + " Max = " + str(max_) + " Min = " + str(min_))
                print(sep_eq_to_loc)

            if max_ > angle_treshold:
                raise Exception("Error: Maximum angle separation = " + str(max_) + " > " + str(angle_treshold))

        if set(['gal_lon','gal_lat']).issubset(table_read.columns):
            
            sep_gal_to_loc = test_benchmark_conversion(table_loc_to_gal, table_gal_to_loc)
            
            mean_ = sep_gal_to_loc.mean()
            min_ = sep_gal_to_loc.min()
            max_ = sep_gal_to_loc.max()
            
            if to_print == True:
                print("GAL TO LOC")            
                print(sep_gal_to_loc)
                print("Mean = " + str(mean_) + " Max = " + str(max_) + " Min = " + str(min_))

            if max_ > angle_treshold:
                raise Exception("Error: Maximum angle separation = " + str(max_) + " > " + str(angle_treshold))

    print("\nChecking Az/Ze")
    if set(['azimuth','zenith']).issubset(table_read.columns):

        if set(['phi','theta']).issubset(table_read.columns):
            sep_loc_to_utm = test_benchmark_conversion(table_utm_to_loc, table_loc_to_utm)
            
            mean_ = sep_loc_to_utm.mean()
            min_ = sep_loc_to_utm.min()
            max_ = sep_loc_to_utm.max()
            
            if to_print == True:
                print("LOC TO UTM")
                print(sep_loc_to_utm)
                print("Mean = " + str(mean_) + " Max = " + str(max_) + " Min = " + str(min_))

            if max_ > angle_treshold:
                raise Exception("Error: Maximum angle separation = " + str(max_) + " > " + str(angle_treshold))


        if set(['RA-J2000','DEC-J2000']).issubset(table_read.columns):

            sep_eq_to_utm = test_benchmark_conversion(table_utm_to_eq, table_eq_to_utm)

            mean_ = sep_eq_to_utm.mean()
            min_ = sep_eq_to_utm.min()
            max_ = sep_eq_to_utm.max()

            if to_print == True:
                print("EQ TO UTM")
                print(sep_eq_to_utm)
                print("Mean = " + str(mean_) + " Max = " + str(max_) + " Min = " + str(min_))
            
            if max_ > angle_treshold:
                raise Exception("Error: Maximum angle separation = " + str(max_) + " > " + str(angle_treshold))


        if set(['gal_lon','gal_lat']).issubset(table_read.columns):
            sep_gal_to_utm = test_benchmark_conversion(table_utm_to_gal, table_gal_to_utm)

            mean_ = sep_gal_to_utm.mean()
            min_ = sep_gal_to_utm.min()
            max_ = sep_gal_to_utm.max()

            if to_print == True:
                print("GAL TO UTM")
                print(sep_gal_to_utm)
                print("Mean = " + str(mean_) + " Max = " + str(max_) + " Min = " + str(min_))

            if max_ > angle_treshold:
                raise Exception("Error: Maximum angle separation = " + str(max_) + " > " + str(angle_treshold))


    print("\nChecking RA/DEC")
    if set(['RA-J2000','DEC-J2000']).issubset(table_read.columns):

        if set(['phi','theta']).issubset(table_read.columns):
            sep_loc_to_eq = test_benchmark_conversion(table_eq_to_loc, table_loc_to_eq)

            mean_ = sep_loc_to_eq.mean()
            min_ = sep_loc_to_eq.min()
            max_ = sep_loc_to_eq.max()

            if to_print == True:
                print("LOC TO EQ")
                print(sep_loc_to_eq)
                print("Mean = " + str(mean_) + " Max = " + str(max_) + " Min = " + str(min_))

            if max_ > angle_treshold:
                raise Exception("Error: Maximum angle separation = " + str(max_) + " > " + str(angle_treshold))


        if set(['azimuth','zenith']).issubset(table_read.columns):
            sep_utm_to_eq = test_benchmark_conversion(table_eq_to_utm, table_utm_to_eq)

            mean_ = sep_utm_to_eq.mean()
            min_ = sep_utm_to_eq.min()
            max_ = sep_utm_to_eq.max()

            if to_print == True:
                print("UTM TO EQ")
                print(sep_utm_to_eq)
                print("Mean = " + str(mean_) + " Max = " + str(max_) + " Min = " + str(min_))

            if max_ > angle_treshold:
                raise Exception("Error: Maximum angle separation = " + str(max_) + " > " + str(angle_treshold))

        if set(['gal_lon','gal_lat']).issubset(table_read.columns):
            sep_gal_to_eq = test_benchmark_conversion(table_eq_to_gal, table_gal_to_eq)

            mean_ = sep_gal_to_eq.mean()
            min_ = sep_gal_to_eq.min()
            max_ = sep_gal_to_eq.max()

            if to_print == True:
                print("GAL TO EQ")
                print(sep_gal_to_eq)
                print("Mean = " + str(mean_) + " Max = " + str(max_) + " Min = " + str(min_))

            if max_ > angle_treshold:
                raise Exception("Error: Maximum angle separation = " + str(max_) + " > " + str(angle_treshold))

    print("\nChecking GAL L/B")
    if set(['gal_lon','gal_lat']).issubset(table_read.columns):

        if set(['phi','theta']).issubset(table_read.columns):
            sep_loc_to_gal = test_benchmark_conversion(table_gal_to_loc, table_loc_to_gal)

            mean_ = sep_loc_to_gal.mean()
            min_ = sep_loc_to_gal.min()
            max_ = sep_loc_to_gal.max()
            
            if to_print == True:
                print("LOC TO GAL")
                print(sep_loc_to_gal)
                print("Mean = " + str(mean_) + " Max = " + str(max_) + " Min = " + str(min_))

            if max_ > angle_treshold:
                raise Exception("Error: Maximum angle separation = " + str(max_) + " > " + str(angle_treshold))

        if set(['azimuth','zenith']).issubset(table_read.columns):
            sep_utm_to_gal = test_benchmark_conversion(table_gal_to_utm, table_utm_to_gal)

            mean_ = sep_utm_to_gal.mean()
            min_ = sep_utm_to_gal.min()
            max_ = sep_utm_to_gal.max()

            if to_print == True:
                print("UTM TO GAL")
                print(sep_utm_to_gal)
                print("Mean = " + str(mean_) + " Max = " + str(max_) + " Min = " + str(min_))

            if max_ > angle_treshold:
                raise Exception("Error: Maximum angle separation = " + str(max_) + " > " + str(angle_treshold))

        if set(['RA-J2000','DEC-J2000']).issubset(table_read.columns):
            sep_eq_to_gal = test_benchmark_conversion(table_gal_to_eq, table_eq_to_gal)

            mean_ = sep_eq_to_gal.mean()
            min_ = sep_eq_to_gal.min()
            max_ = sep_eq_to_gal.max()

            if to_print == True:
                print("EQ TO GAL")
                print(sep_eq_to_gal)
                print("Mean = " + str(mean_) + " Max = " + str(max_) + " Min = " + str(min_))

            if max_ > angle_treshold:
                raise Exception("Error: Maximum angle separation = " + str(max_) + " > " + str(angle_treshold))

    print("... End of Angle Separation Test")
    return 0

def test_separation(file0 = "", detector_ = "antares", detector_to_ = "antares"):

    print("Starting Calculated value Separation Test...")
    if file0 == "":
        file0 = "/home/htedjditi/work/test_astro/antares_coordinate_systems_benchmark.csv"

    print("Read " + file0 + " file. \n")

    table_read = reader_from_file(file0)
    print(table_read)

    atol_value = 0.021

    if set(['RA-J2000','DEC-J2000']).issubset(table_read.columns):
        true_ra = table_read['RA-J2000']
        if type(true_ra[0]) == str:
            true_ra = Angle(true_ra, unit = 'hourangle')
        elif type(true_ra[0]) == np.float64:
            true_ra = Angle(true_ra, unit = u.deg)

        true_dec = table_read['DEC-J2000']
        if type(true_dec[0]) == str:
            true_dec = Angle(true_dec, unit = u.deg)
        elif type(true_dec[0]) == np.float64:
            true_dec = Angle(true_dec, unit = u.deg)
        
    print("Asserting from Phi, theta")
    if set(['phi','theta']).issubset(table_read.columns):

        table_loc_to_utm = transform_to_new_frame(table_read, "Detector", "UTM", detector_, detector_to_)
        table_loc_to_eq = transform_to_new_frame(table_read, "Detector", "equatorial", detector_, detector_to_)
        table_loc_to_gal = transform_to_new_frame(table_read, "Detector", "galactic", detector_, detector_to_)

        phi_, theta_ = zip(*table_loc_to_utm.apply( lambda x: get_phi_theta(x.SkyCoord_base, detector_to_) , axis = 1))

        az_, ze_ = zip(*table_loc_to_utm.apply( lambda x: get_az_zenith(x.SkyCoord_new, detector_to_) , axis = 1))
        
        ra_, dec_ = zip(*table_loc_to_eq.apply( lambda x: get_ra_dec(x.SkyCoord_new) , axis = 1))     
        l_, b_ = zip(*table_loc_to_gal.apply( lambda x: get_l_b(x.SkyCoord_new) , axis = 1))     

        if set(['phi','theta']).issubset(table_read.columns):
            np.testing.assert_allclose(table_read['phi'],phi_,rtol=0,atol=atol_value,verbose=True)
            np.testing.assert_allclose(table_read['theta'],theta_,rtol=0,atol=atol_value,verbose=True)

        if set(['azimuth','zenith']).issubset(table_read.columns):
            np.testing.assert_allclose(table_read['azimuth'],az_,rtol=0,atol=atol_value,verbose=True)
            np.testing.assert_allclose(table_read['zenith'],ze_,rtol=0,atol=atol_value,verbose=True)
        
        if set(['RA-J2000','DEC-J2000']).issubset(table_read.columns):
            np.testing.assert_allclose(true_ra.deg,ra_,rtol=0,atol=atol_value,verbose=True)
            np.testing.assert_allclose(true_dec.deg,dec_,rtol=0,atol=atol_value,verbose=True)
        
        if set(['gal_lon','gal_lat']).issubset(table_read.columns):
            np.testing.assert_allclose(table_read['gal_lon'],l_,rtol=0,atol=atol_value,verbose=True)
            np.testing.assert_allclose(table_read['gal_lat'],b_,rtol=0,atol=atol_value,verbose=True)

    print("Asserting from Az, Zenith")
    if set(['azimuth','zenith']).issubset(table_read.columns):

        table_utm_to_loc = transform_to_new_frame(table_read, "UTM", "Detector", detector_, detector_to_)
        table_utm_to_eq = transform_to_new_frame(table_read, "UTM", "equatorial", detector_, detector_to_)
        table_utm_to_gal = transform_to_new_frame(table_read, "UTM", "galactic", detector_, detector_to_)

        az_, ze_ = zip(*table_utm_to_loc.apply( lambda x: get_az_zenith(x.SkyCoord_base, detector_to_) , axis = 1))

        phi_, theta_ = zip(*table_utm_to_loc.apply( lambda x: get_phi_theta(x.SkyCoord_new, detector_to_) , axis = 1))
        
        ra_, dec_ = zip(*table_utm_to_eq.apply( lambda x: get_ra_dec(x.SkyCoord_new) , axis = 1))     
        l_, b_ = zip(*table_utm_to_gal.apply( lambda x: get_l_b(x.SkyCoord_new) , axis = 1))     

        if set(['phi','theta']).issubset(table_read.columns):
            np.testing.assert_allclose(table_read['phi'],phi_,rtol=0,atol=atol_value,verbose=True)
            np.testing.assert_allclose(table_read['theta'],theta_,rtol=0,atol=atol_value,verbose=True)

        if set(['azimuth','zenith']).issubset(table_read.columns):
            np.testing.assert_allclose(table_read['azimuth'],az_,rtol=0,atol=atol_value,verbose=True)
            np.testing.assert_allclose(table_read['zenith'],ze_,rtol=0,atol=atol_value,verbose=True)
        
        if set(['RA-J2000','DEC-J2000']).issubset(table_read.columns):
            np.testing.assert_allclose(true_ra.deg,ra_,rtol=0,atol=atol_value,verbose=True)
            np.testing.assert_allclose(true_dec.deg,dec_,rtol=0,atol=atol_value,verbose=True)
        
        if set(['gal_lon','gal_lat']).issubset(table_read.columns):
            np.testing.assert_allclose(table_read['gal_lon'],l_,rtol=0,atol=atol_value,verbose=True)
            np.testing.assert_allclose(table_read['gal_lat'],b_,rtol=0,atol=atol_value,verbose=True)

    print("Asserting from ra, dec")
    if set(['RA-J2000','DEC-J2000']).issubset(table_read.columns):

        table_eq_to_utm = transform_to_new_frame(table_read, "equatorial", "UTM", detector_, detector_to_)
        table_eq_to_loc = transform_to_new_frame(table_read, "equatorial", "Detector", detector_, detector_to_)
        table_eq_to_gal = transform_to_new_frame(table_read, "equatorial", "galactic", detector_, detector_to_)

        ra_, dec_ = zip(*table_eq_to_utm.apply( lambda x: get_ra_dec(x.SkyCoord_base) , axis = 1))     
        
        phi_, theta_ = zip(*table_eq_to_loc.apply( lambda x: get_phi_theta(x.SkyCoord_new, detector_to_) , axis = 1))
        az_, ze_ = zip(*table_eq_to_utm.apply( lambda x: get_az_zenith(x.SkyCoord_new, detector_to_) , axis = 1))
        l_, b_ = zip(*table_eq_to_gal.apply( lambda x: get_l_b(x.SkyCoord_new) , axis = 1))     

        if set(['phi','theta']).issubset(table_read.columns):
            np.testing.assert_allclose(table_read['phi'],phi_,rtol=0,atol=atol_value,verbose=True)
            np.testing.assert_allclose(table_read['theta'],theta_,rtol=0,atol=atol_value,verbose=True)

        if set(['azimuth','zenith']).issubset(table_read.columns):
            np.testing.assert_allclose(table_read['azimuth'],az_,rtol=0,atol=atol_value,verbose=True)
            np.testing.assert_allclose(table_read['zenith'],ze_,rtol=0,atol=atol_value,verbose=True)
        
        if set(['RA-J2000','DEC-J2000']).issubset(table_read.columns):
            np.testing.assert_allclose(true_ra.deg,ra_,rtol=0,atol=atol_value,verbose=True)
            np.testing.assert_allclose(true_dec.deg,dec_,rtol=0,atol=atol_value,verbose=True)
        
        if set(['gal_lon','gal_lat']).issubset(table_read.columns):
            np.testing.assert_allclose(table_read['gal_lon'],l_,rtol=0,atol=atol_value,verbose=True)
            np.testing.assert_allclose(table_read['gal_lat'],b_,rtol=0,atol=atol_value,verbose=True)

    print("Asserting from l, b")
    if set(['gal_lon','gal_lat']).issubset(table_read.columns):

        table_gal_to_loc = transform_to_new_frame(table_read, "galactic", "Detector", detector_, detector_to_)
        table_gal_to_eq = transform_to_new_frame(table_read, "galactic", "equatorial", detector_, detector_to_)
        table_gal_to_utm = transform_to_new_frame(table_read, "galactic", "UTM", detector_, detector_to_)

        l_, b_ = zip(*table_gal_to_utm.apply( lambda x: get_l_b(x.SkyCoord_base) , axis = 1))     
                
        phi_, theta_ = zip(*table_utm_to_loc.apply( lambda x: get_phi_theta(x.SkyCoord_new, detector_to_) , axis = 1))
        az_, ze_ = zip(*table_eq_to_utm.apply( lambda x: get_az_zenith(x.SkyCoord_new, detector_to_) , axis = 1))

        ra_, dec_ = zip(*table_utm_to_eq.apply( lambda x: get_ra_dec(x.SkyCoord_new) , axis = 1))

        if set(['phi','theta']).issubset(table_read.columns):
            np.testing.assert_allclose(table_read['phi'],phi_,rtol=0,atol=atol_value,verbose=True)
            np.testing.assert_allclose(table_read['theta'],theta_,rtol=0,atol=atol_value,verbose=True)

        if set(['azimuth','zenith']).issubset(table_read.columns):
            np.testing.assert_allclose(table_read['azimuth'],az_,rtol=0,atol=atol_value,verbose=True)
            np.testing.assert_allclose(table_read['zenith'],ze_,rtol=0,atol=atol_value,verbose=True)
        
        if set(['RA-J2000','DEC-J2000']).issubset(table_read.columns):
            np.testing.assert_allclose(true_ra.deg,ra_,rtol=0,atol=atol_value,verbose=True)
            np.testing.assert_allclose(true_dec.deg,dec_,rtol=0,atol=atol_value,verbose=True)
        
        if set(['gal_lon','gal_lat']).issubset(table_read.columns):
            np.testing.assert_allclose(table_read['gal_lon'],l_,rtol=0,atol=atol_value,verbose=True)
            np.testing.assert_allclose(table_read['gal_lat'],b_,rtol=0,atol=atol_value,verbose=True)


    print("... End of Calculated Value Separation Test")
    return 0

def main():
    #test_coordinate_transformator()
    #test_coordinate_transformator(data_path("astro/ARCA_coordinate_systems_benchmark.csv"), "arca", "arca")
    #test_coordinate_transformator(data_path("astro/ORCA_coordinate_systems_benchmark.csv"), "orca", "orca")
    #test_coordinate_transformator(data_path("astro/ARCA_astro_objects_benchmark.csv"), "arca", "arca", type_name = "astro_objects")
    #test_coordinate_transformator(data_path("astro/ORCA_moon_sun_position_benchmark.csv"), "orca", "orca", type_name = "moon_sun_position")

    #testing_angle_separation
    
    #coordinate system benchmark
    #ok#
    test_angle_separation(to_print = True)
    #ok#test_angle_separation(data_path("astro/ARCA_coordinate_systems_benchmark.csv"), "arca", "arca")
    #ok#test_angle_separation(data_path("astro/ORCA_coordinate_systems_benchmark.csv"), "orca", "orca")

    #astro object benchmark
    #ok#test_angle_separation(data_path("astro/antares_astro_objects_benchmark.csv"), "antares", "antares")
    #ok#test_angle_separation(data_path("astro/ARCA_astro_objects_benchmark.csv"), "arca", "arca")
    #ok#test_angle_separation(data_path("astro/ORCA_astro_objects_benchmark.csv"), "orca", "orca")
    
    
    #moon sun postion benchmark
    #ok#test_angle_separation(data_path("astro/antares_moon_sun_position_benchmark.csv"), "antares", "antares")
    #no#
    #test_angle_separation(data_path("astro/ARCA_moon_sun_position_benchmark.csv"), "arca", "arca")
    #no#test_angle_separation(data_path("astro/ORCA_moon_sun_position_benchmark.csv"), "orca", "orca")
    

    #test_separation value per value
    
    #Working#

    #astro_object
    #test_separation(data_path("astro/antares_astro_objects_benchmark.csv"), "antares", "antares")
    #test_separation(data_path("astro/ARCA_astro_objects_benchmark.csv"), "arca", "arca")
    #test_separation(data_path("astro/ORCA_astro_objects_benchmark.csv"), "orca", "orca")

    #moon_sun_position
    #test_separation(data_path("astro/antares_moon_sun_position_benchmark.csv"), "antares", "antares")

    #NotWorking#
    
    #coordinate benchmark

    #2 degree max
    #test_separation()
    #0.9 degree max"
    #test_separation(data_path("astro/ARCA_coordinate_systems_benchmark.csv"), "arca", "arca")
    #0.95 degree max
    #test_separation(data_path("astro/ORCA_coordinate_systems_benchmark.csv"), "orca", "orca")

    #sun moon position
    #totally wrong for moon
    #test_separation(data_path("astro/ARCA_moon_sun_position_benchmark.csv"), "arca", "arca")
    #totally wrong for moon
    #test_separation(data_path("astro/ORCA_moon_sun_position_benchmark.csv"), "orca", "orca")

    
if __name__ == "__main__":
    main()
