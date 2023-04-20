import os
import platform
import re
import subprocess
import sys
from pathlib import Path

from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext

with open("splinepy/_version.py") as f:
    version = eval(f.read().strip().split("=")[-1])

# taken from pybind cmake example
# https://github.com/pybind/cmake_example/blob/master/setup.py

# Convert distutils Windows platform specifiers to CMake -A arguments
PLAT_TO_CMAKE = {
    "win32": "Win32",
    "win-amd64": "x64",
    "win-arm32": "ARM",
    "win-arm64": "ARM64",
}


# A CMakeExtension needs a sourcedir instead of a file list.
# The name must be the _single_ output extension from the CMake build.
# If you need multiple extensions, see scikit-build.
class CMakeExtension(Extension):
    def __init__(
        self, name: str, sourcedir: str = "", extra_args: dict = None
    ) -> None:
        super().__init__(name, sources=[])
        self.sourcedir = os.fspath(Path(sourcedir).resolve())
        if extra_args is not None:
            self.extra_args = extra_args
            self.debug = extra_args.get("--debug", False)


class CMakeBuild(build_ext):
    def build_extension(self, ext: CMakeExtension) -> None:
        # Must be in this form due to bug in .resolve()
        # only fixed in Python 3.10+
        ext_fullpath = Path.cwd() / self.get_ext_fullpath(
            ext.name
        )  # type: ignore[no-untyped-call]
        extdir = ext_fullpath.parent.resolve()

        # Using this requires trailing slash for auto-detection & inclusion of
        # auxiliary "native" libs

        debug = (
            int(os.environ.get("DEBUG", 0))
            if self.debug is None
            else self.debug
        )
        # overwrite if ext.debug exists
        debug = ext.debug if hasattr(ext, "debug") else debug
        cfg = "Debug" if debug else "Release"

        # CMake lets you override the generator - we need to check this.
        # Can be set with Conda-Build, for example.
        cmake_generator = os.environ.get("CMAKE_GENERATOR", "")

        # Set Python_EXECUTABLE instead if you use PYBIND11_FINDPYTHON
        # EXAMPLE_VERSION_INFO shows you how to pass a value into the C++ code
        # from Python.
        cmake_args = [
            f"-DCMAKE_LIBRARY_OUTPUT_DIRECTORY={extdir}{os.sep}",
            f"-DPYTHON_EXECUTABLE={sys.executable}",
            f"-DCMAKE_BUILD_TYPE={cfg}",  # not used on MSVC, but no harm
        ]
        # extra cmake args
        cmake_args.extend(ext.extra_args["cmake_args"])

        build_args = []
        # Adding CMake arguments set as environment variable
        # (needed e.g. to build for ARM OSx on conda-forge)
        if "CMAKE_ARGS" in os.environ:
            cmake_args += [
                item for item in os.environ["CMAKE_ARGS"].split(" ") if item
            ]

        if self.compiler.compiler_type != "msvc":
            # Using Ninja-build since it a) is available as a wheel and b)
            # multithreads automatically. MSVC would require all variables be
            # exported for Ninja to pick it up, which is a little tricky to do.
            # Users can override the generator with CMAKE_GENERATOR in CMake
            # 3.15+.
            if not cmake_generator or cmake_generator == "Ninja":
                try:
                    import ninja  # noqa: F401

                    ninja_executable_path = Path(ninja.BIN_DIR) / "ninja"
                    cmake_args += [
                        "-GNinja",
                        "-DCMAKE_MAKE_PROGRAM"
                        f":FILEPATH={ninja_executable_path}",
                    ]
                except ImportError:
                    pass

        else:
            # Single config generators are handled "normally"
            single_config = any(
                x in cmake_generator for x in {"NMake", "Ninja"}
            )

            # CMake allows an arch-in-generator style
            # for backward compatibility
            contains_arch = any(x in cmake_generator for x in {"ARM", "Win64"})

            # Specify the arch if using MSVC generator, but only if it doesn't
            # contain a backward-compatibility arch spec already in the
            # generator name.
            if not single_config and not contains_arch:
                cmake_args += ["-A", PLAT_TO_CMAKE[self.plat_name]]

            # Multi-config generators have a different way to specify configs
            if not single_config:
                cmake_args += [
                    "-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_"
                    f"{cfg.upper()}={extdir}"
                ]
                build_args += ["--config", cfg]

        if sys.platform.startswith("darwin"):
            # Cross-compile support for macOS - respect ARCHFLAGS if set
            archs = re.findall(r"-arch (\S+)", os.environ.get("ARCHFLAGS", ""))
            if archs:
                cmake_args += [
                    "-DCMAKE_OSX_ARCHITECTURES={}".format(";".join(archs))
                ]

        # Set CMAKE_BUILD_PARALLEL_LEVEL to control the parallel build level
        # across all generators.
        if "CMAKE_BUILD_PARALLEL_LEVEL" not in os.environ:
            # self.parallel is a Python 3 only way to set parallel jobs by hand
            # using -j in the build_ext call,
            # not supported by pip or PyPA-build.
            if hasattr(self, "parallel") and self.parallel:
                # CMake 3.12+ only.
                build_args += [f"-j{self.parallel}"]

            if len(ext.extra_args["build_args"]) != 0:
                build_args.extend(ext.extra_args["build_args"])

        build_temp = Path(self.build_temp) / ext.name
        if not build_temp.exists():
            build_temp.mkdir(parents=True)

        subprocess.run(
            ["cmake", ext.sourcedir] + cmake_args, cwd=build_temp, check=True
        )
        subprocess.run(
            ["cmake", "--build", "."] + build_args, cwd=build_temp, check=True
        )


with open("README.md") as readme:
    long_description = readme.read()


def if_in_argv_do(argv, flag, build_opts, value):
    if flag in argv:
        argv.remove(flag)
        if flag in build_opts["cmake_args"]:
            build_opts["cmake_args"][flag] = value
        elif flag in build_opts["build_args"]:
            build_opts["build_args"][flag] = value
        elif flag in build_opts:  # debug
            build_opts[flag] = value
        else:
            raise ValueError(f"Invalid build option flag ({flag})")
        return True
    else:
        return False


# keys / options
keys = (
    "--verbose_make",
    "--minimal",
    "--enable_warning",
    "--serial_build",
    "--debug",
    "--no_explicit",
)

# build options with default options
build_options = {
    "cmake_args": {
        keys[0]: "-DSPLINEPY_VERBOSE_MAKE=OFF",
        keys[1]: "-DSPLINEPY_MORE=ON",
        keys[2]: "-DSPLINEPY_ENABLE_WARNINGS=OFF",
        keys[5]: "-DSPLINEPY_BUILD_EXPLICIT=ON",
    },
    "build_args": {keys[3]: f"-j {os.cpu_count()}"},
    keys[4]: False,
}

# go through options
if_in_argv_do(sys.argv, keys[0], build_options, "-DSPLINEPY_VERBOSE_MAKE=ON")
if_in_argv_do(sys.argv, keys[1], build_options, "-DSPLINEPY_MORE=OFF")
if_in_argv_do(
    sys.argv, keys[2], build_options, "-DSPLINEPY_ENABLE_WARNINGS=ON"
)
if_in_argv_do(sys.argv, keys[3], build_options, "-j 1")
if_in_argv_do(sys.argv, keys[4], build_options, True)
if_in_argv_do(
    sys.argv, keys[5], build_options, "-DSPLINEPY_BUILD_EXPLICIT=OFF"
)

# environ option
if eval(os.environ.get("SPLINEPY_MINIMAL_DEBUG_BUILD", "False")):
    print("Environment variable SPLINEPY_MINIMAL_DEBUG_BUILD set.")
    build_options["cmake_args"][keys[1]] = "-DSPLINEPY_MORE=OFF"
    build_options[keys[4]] = True

if (
    "SPLINEPY_GITHUB_ACTIONS_BUILD" in os.environ
    and platform.system().startswith("Windows")
):
    print("Environment variable SPLINEPY_GITHUB_ACTIONS_BUILD set.")
    print("Platform is detected as Windows.")
    build_options["cmake_args"][keys[1]] = "-DSPLINEPY_MORE=OFF"

# print options and pack items as list
print("************* build options ************")
print("{")
for key, value in build_options.items():
    # print
    if isinstance(value, dict):
        print(f"  {key}: " "{")
        for kk, vv in value.items():
            print(f"    {kk}: {vv}")
        print("  }")
    else:
        print(f"  {key} : {value}")
    # pack as list
    if key.startswith("--debug"):
        continue
    build_options[key] = list(value.values())
print("}")
print("****************************************")


setup(
    name="splinepy",
    version=version,
    author="Jaewook Lee",
    author_email="jaewooklee042@gmail.com",
    description="Python N-Dimensional Bezier, RationalBezier, "
    "BSpline and NURBS library with C++ Backend.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/tataratat/splinepy",
    packages=[
        "splinepy",
        "splinepy.io",
        "splinepy.utils",
        "splinepy.helpme",
    ],
    install_requires=[
        "numpy",
    ],
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Natural Language :: English",
        "Topic :: Scientific/Engineering",
    ],
    ext_modules=[
        CMakeExtension("splinepy.splinepy_core", extra_args=build_options)
    ],
    cmdclass=dict(build_ext=CMakeBuild),
    extras_require={"test": ["pytest>=6.0"]},
    zip_safe=False,
    license="MIT",
)
