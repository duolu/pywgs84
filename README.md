# pywgs84

The **pywgs84** library is a Python library for working with coordinates in WGS84. This is mainly developed for a project that needs to track vehicles using cameras while generating trajectories in latitude and longitude to compare with GPS results. It is also designed to teach and learn basic concepts in WGS84.

It has the following modules.

* [pywgs84.py](https://github.com/duolu/pywgs84/blob/master/pywgs84.py) - The core functions for working with coordinates in WGS84.

The code is tested on Python 3.6 with NumPy 1.19.

## Terminology

* "**llh**" or "**LLH**" means **(latitude, longitude, height)**.
* "latitude" and "longitude" are always in degrees.
* "lat" and "lon" are always in radians (only used internally).
* "height" is the height above the geodetic ellipsoid, in meters.
* "llh_ref" is a tuple of LLH for specifying a reference point.

* "**ecef**" or "**ECEF**" means "**Earth-Centered-Earth-Fixed**".
* "**xyz**" or "**XYZ**" means coordinates in the **ECEF** reference frame.
* "**enu**" or "**ENU**" means **(east, north, up)**.
* "**ned**" or "**NED**" means **(north, east, down)**.
* "**esd**" or "**ESD**" means **(east, south, down)**.

NOTE: All Euclidean reference frames are right handed. All Euclidean coordinates are in meters. "height" is not "altitude". "height" is the displacement relative to the **geodetic ellipsoid** in WGS84, not the **geoid** or the actual surface of the earth. There is a reference point needed for ENU/NED/ESD frame, i.e., the point on the ellipsoid to set up a tangent plane.


## Currently Implemented Functions

The following functions are implemented currently.

* Conversions between coordinates in LLH and ECEF.
* Conversions between coordinates in ECEF and ENU/NED/ESD. Hence, it can also do conversions between LLH and ENU/NED/ESD.

There is not much code and it essentially explains itself. More features on geodesy will be added as needded.

