import numpy as np

from astropy.coordinates import ICRS, Galactic, FK4, FK5  # Low-level frames
import astropy.units as u
from astropy.coordinates import Angle

from astropy.coordinates import BaseCoordinateFrame
import astropy.coordinates.representation as rep
from astropy.coordinates import RepresentationMapping
from astropy.coordinates.attributes import TimeAttribute, EarthLocationAttribute, QuantityAttribute
from astropy.coordinates import TransformGraph, frame_transform_graph, FunctionTransform

from km3astro.coord import convergence_angle


class Detector(BaseCoordinateFrame):
    default_representation = rep.PhysicsSphericalRepresentation

    #Specify frame attributes required to fully specify the frame
    obstime = TimeAttribute(default=None)
    location = EarthLocationAttribute(default=get_location("arca"))


class UTM(BaseCoordinateFrame):

    default_representation = rep.PhysicsSphericalRepresentation
    
    frame_specific_representation_info = {
        rep.PhysicsSphericalRepresentation: [RepresentationMapping('phi', 'azimuth'),
                                    RepresentationMapping('theta', 'zenith')]
    }

    obstime = TimeAttribute(default=None)
    location = EarthLocationAttribute(default=get_location("arca"))
    
@frame_transform_graph.transform(FunctionTransform, Detector, AltAz)
def Detector_to_AltAz(Detector_, altaz):

    phi = Detector_.phi.rad
    theta = Detector_.theta.rad
    loc = Detector_.location
    time = Detector_.obstime

    conv_angle = Angle(convergence_angle(loc.lat.rad, loc.lon.rad), unit = u.radian)

    altitude = theta - np.pi / 2
    Corrected_azimuth = (np.pi / 2 - phi + np.pi + conv_angle.rad) % (2 * np.pi)

    altaz = AltAz(alt = altitude *rad, az = Corrected_azimuth *rad, obstime=time, location=loc)

    return altaz

@frame_transform_graph.transform(FunctionTransform, AltAz, Detector)
def AltAz_to_Detector(altaz_, detector_):

    alt = altaz_.alt.rad
    az = altaz_.az.rad
    loc = altaz_.location
    time = altaz_.obstime
    r = u.Quantity(100, u.m) # Warning dummy r value! For SphericalRepresentation to PhysicsSphericalRepresentation conversion
    
    conv_angle = Angle(convergence_angle(loc.lat.rad, loc.lon.rad), unit = u.radian)

    phi = np.pi/2 + np.pi + conv_angle.rad - az
    
    theta = alt + np.pi/2

    detector_ = Detector(phi = phi *rad, theta = theta *rad, obstime=time, location=loc, r = r)

    return detector_


@frame_transform_graph.transform(FunctionTransform, UTM, AltAz)
def UTM_to_AltAz(UTM_, altaz):

    az = UTM_.azimuth.rad
    ze = UTM_.zenith.rad

    loc = UTM_.location
    time = UTM_.obstime

    conv_angle = Angle(convergence_angle(loc.lat.rad, loc.lon.rad), unit = u.radian)

    altitude = np.pi / 2 - ze
    Corrected_azimuth = (np.pi / 2 - az + 2*np.pi + conv_angle.rad) % (2 * np.pi)

    altaz = AltAz(alt = altitude *rad, az = Corrected_azimuth *rad, obstime=time, location=loc)

    return altaz

@frame_transform_graph.transform(FunctionTransform, AltAz, UTM)
def AltAz_to_UTM(altaz_, UTM_):

    alt = altaz_.alt.rad
    caz = altaz_.az.rad

    loc = altaz_.location
    time = altaz_.obstime
    r = u.Quantity(100, u.m) # Warning dummy r value! For SphericalRepresentation to PhysicsSphericalRepresentation.
    
    conv_angle = Angle(convergence_angle(loc.lat.rad, loc.lon.rad), unit = u.radian)

    az = 5*np.pi/2 + conv_angle.rad - caz
    az = az %(2*np.pi)
    ze = np.pi/2 - alt

    UTM_ = UTM(azimuth = az *rad, zenith = ze *rad, obstime = time, location = loc, r = r)
    return UTM_

@frame_transform_graph.transform(FunctionTransform, UTM, Detector)
def UTM_to_Detector(UTM_, detector):
    
    phi = UTM_.azimuth.rad - np.pi
    theta = np.pi - UTM_.zenith.rad
    loc = UTM_.location
    time = UTM_.obstime
    r = UTM_.r

    detector = Detector(phi = phi *rad, theta = theta *rad, obstime=time, location=loc, r = r)

    return detector

@frame_transform_graph.transform(FunctionTransform, Detector, UTM)
def Detector_to_UTM(detector, UTM_):
    
    az = detector.phi.rad + np.pi
    ze = np.pi - detector.theta.rad

    loc = detector.location
    time = detector.obstime
    r = detector.r

    UTM_ = UTM(azimuth = az *rad, zenith = ze *rad, obstime = time, location = loc, r = r)

    return UTM_
