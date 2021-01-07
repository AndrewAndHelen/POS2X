# pos2ecf

a tool which could transfer POS to ECEF


## Dependencies

* Eigen

##  Compilation
prerequisites: cmake version>=3.0

```
1. git clone https://github.com/AndrewAndJenny/pos2ecf.git
2. cd pos2ecf
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
* The format of image
**(&lambda;, &omega;, H, heading, pitch, roll)**
>*Tip：(&lambda;, &omega;, H) is object locations in lat–long–height.*  

>*Tip：(heading, pitch, roll) is orientations in a heading–pitch–roll format follow NED coordinate system*

## Contact

If you found bugs or have new ideas,please pull requests😄  

If you have trouble compiling or using this software,email to [15313326374@163.com](mailto:15313326374@163.com)

## Good Luck For You!
