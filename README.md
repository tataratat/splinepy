# splinepy
Python BSpline and library with C++ backend, [SplineLib](https://github.com/SplineLib/SplineLib).


## Install-guide
Option 1: `pip`
```
pip install splinepy
```

Option 2: `setup.py`  
You can locally compile and install using `setup.py`.
This requires C++20 (C++17 for release mode) compatible C++ compiler
and cmake version 3.16 or higher.
```
git clone git@github.com:tataratat/splinepy.git
cd splinepy
git submodule update --init
python3 setup.py install
```
