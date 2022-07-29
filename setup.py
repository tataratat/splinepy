import os
import re
import sys
import platform
import subprocess

from distutils.version import LooseVersion
from setuptools import setup, find_packages, Extension
from setuptools.command.build_ext import build_ext


with open("splinepy/_version.py") as f:
    version = eval(f.read().strip().split("=")[-1])


class CMakeExtension(Extension):
    def __init__(self, name, sourcedir='', cmake_args=None):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)
        if cmake_args is not None:
            self.cmake_args = cmake_args


class CMakeBuild(build_ext):
    def run(self):
        try:
            out = subprocess.check_output(['cmake', '--version'])
        except OSError:
            raise RuntimeError(
                "CMake must be installed to build the following extensions: " +
                ", ".join(e.name for e in self.extensions))

        if platform.system() == "Windows":
            cmake_version = LooseVersion(re.search(r'version\s*([\d.]+)',
                                         out.decode()).group(1))
            if cmake_version < '3.1.0':
                raise RuntimeError("CMake >= 3.1.0 is required on Windows")

        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):
        extdir = os.path.abspath(
            os.path.dirname(self.get_ext_fullpath(ext.name)))
        cmake_args = ['-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=' + extdir,
                      '-DPYTHON_EXECUTABLE=' + sys.executable]

        # runtime cmake args
        cmake_args += ext.cmake_args

        cfg = 'Debug' if self.debug else 'Release'
        build_args = ['--config', cfg]

        if platform.system() == "Windows":
            cmake_args += ['-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{}={}'.format(
                cfg.upper(),
                extdir)]
            if sys.maxsize > 2**32:
                cmake_args += ['-A', 'x64']
            build_args += ['--', '/m']
        else:
            cmake_args += ['-DCMAKE_BUILD_TYPE=' + cfg]
            build_args += ['--', '-j4']

        env = os.environ.copy()
        env['CXXFLAGS'] = '{} -DVERSION_INFO=\\"{}\\"'.format(
            env.get('CXXFLAGS', ''),
            self.distribution.get_version())
        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)
        subprocess.check_call(['cmake', ext.sourcedir] + cmake_args,
                              cwd=self.build_temp, env=env)
        subprocess.check_call(['cmake', '--build', '.'] + build_args,
                              cwd=self.build_temp)
        print()  # Add an empty line for cleaner output


with open("README.md", "r") as readme:
    long_description = readme.read()


# cmake args
# TODO: prettier solution
flags = dict(
    have_splinelib="--have_splinelib",
    verbose_make="--verbose_make",
    minimal="--minimal",
    enable_warning="--enable_warning",
)
cma = []

if flags["have_splinelib"] not in sys.argv:
    print("*** compiling splinelib ***")
    cma.append("-DSPLINEPY_COMPILE_SPLINELIB=ON")
else:
    print("*** NOT compiling splinelib ***")
    sys.argv.remove(flags["have_splinelib"])

if flags["verbose_make"] in sys.argv:
    print("*** verbose make ***")
    cma.append("-DVERBOSE_MAKE=ON")
    sys.argv.remove(flags["verbose_make"])

if flags["minimal"] in sys.argv:
    print("*** compiling only a minimal set of splines")
    cma += ["-DMINIMAL=ON"]
    sys.argv.remove(flags["minimal"])

if flags["enable_warning"] in sys.argv:
    print("*** adding warning flags ***")
    cma += ["-DENABLE_WARNINGS=ON"]
    sys.argv.remove(flags["enable_warning"])


setup(
    name='splinepy',
    version=version,
    author='Jaewook Lee',
    author_email='jlee@ilsb.tuwien.ac.at',
    description='Python NURBS & BSpline library with C++ Backend.',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/tataratat/splinepy",
    packages=[
        "splinepy",
        "splinepy.io",
    ],
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Natural Language :: English',
        'Topic :: Scientific/Engineering',
    ],
    ext_modules=[CMakeExtension('splinepy._splinepy', cmake_args=cma)],
    cmdclass=dict(build_ext=CMakeBuild),
    zip_safe=False,
    license="MIT",
)
