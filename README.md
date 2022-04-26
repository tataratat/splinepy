# splinepy
Python BSpline and NURBS library with C++ backend, [SplineLib](https://github.com/SplineLib/SplineLib).


## Install-guide
Option 1: `pip`
```
pip install --upgrade pip
pip install splinepy
```

Option 2: `setup.py`  
You can locally compile and install using `setup.py`.
This requires C++20 (C++17 for release mode) compatible C++ compiler
and cmake version 3.16 or higher.
```
git clone git@github.com:tataratat/splinepy.git
cd splinepy
python3 setup.py install
```

if there are problems with the vtk version not finding lib9openh264.so.5 create a soft link in the anaconda lib
directory
```
cd <anaconda-directory>/envs/<your conda env>/lib/
ln -s libopenh264.so.2.1.1 libopenh264.so.5
```
