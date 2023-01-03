# POS2X

a tool which could transfer POS


## Dependencies

* Eigen

##  Compilation
prerequisites: cmake version>=3.0

```
1. git clone https://github.com/AndrewAndHelen/POS2X.git
2. cd POS2X
```

**in the cmake-bash**
* Windows
```
1. mkdir build && cd build
2. cmake .. -A x64
3. open the vs solution, and install
```
* Linux
```
1. mkdir build && cd build
2. cmake ..
3. make
```

## Functions
- [x] convert WGS84 longitude, latitude, ellipsoidal height, NED heading pitch roll to WGS84 X, Y, Z, ECEF heading pitch roll
- [x] convert WGS84 longitude, latitude, ellipsoidal height, NED heading pitch roll to WGS84 X, Y, Z, ECEF yaw, pitch, roll
- [x] convert WGS84 longitude, latitude, ellipsoidal height, NED yaw, pitch, roll to CGCS2000(GAUSS) X, Y, Z, PHOTOGRAMMETRY  omega, phi, kappa
- [x] convert WGS84 longitude, latitude, ellipsoidal height, NED yaw, pitch, roll to CGCS2000(UTM) X, Y, Z, PHOTOGRAMMETRY  omega, phi, kappa
- [x] convert WGS84 omega, phi, kappa to ECEF omega, phi, kappa
- [x] convert SMART3D omega, phi, kappa to PHOTOGRAMMETRY omega, phi, kappa

## Contact

If you found bugs or have new ideas,please pull requests😄  

If you have trouble compiling or using this software,email to [15313326374@163.com](mailto:15313326374@163.com)

## Good Luck For You!
