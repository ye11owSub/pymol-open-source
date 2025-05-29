#!/usr/bin/env python
#
# This script only applies if you are performing a Python setuptools-based
# installation of PyMOL.
#
# It may assume that all of PyMOL's external dependencies are
# pre-installed into the system.

import argparse
import glob
import os
import re
import shutil
import sys
import sysconfig
import time
from collections import defaultdict
from itertools import chain
from pathlib import Path
from subprocess import PIPE, Popen

from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext
from setuptools.command.build_py import build_py
from setuptools.command.install import install

# non-empty DEBUG variable turns off optimization and adds -g flag
DEBUG = bool(os.getenv("DEBUG", False))
WIN = sys.platform.startswith("win")
MAC = sys.platform.startswith("darwin")


# Have to copy from "create_shadertext.py" script due to the use of pyproject.toml
# Full explanation:
# https://github.com/pypa/setuptools/issues/3939
def create_all(generated_dir: str, pymol_dir: str = ".") -> None:
    """
    Generate various stuff
    """
    generated_dir_path = Path(generated_dir)
    pymol_dir_path = Path(pymol_dir)

    generated_dir_path.mkdir(parents=True, exist_ok=True)
    pymol_dir_path.mkdir(parents=True, exist_ok=True)

    create_shadertext(
        shader_dir=generated_dir_path / "data" / "shaders",
        shader_dir2=pymol_dir_path,
        output_header=generated_dir_path / "ShaderText.h",
        output_source=generated_dir_path / "ShaderText.cpp",
    )
    create_buildinfo(generated_dir, pymol_dir)


def create_shadertext(
    shader_dir: Path,
    shader_dir2: Path,
    output_header: Path,
    output_source: Path,
) -> None:
    varname = "_shader_cache_raw"
    include_deps = defaultdict(set)
    ifdef_deps = defaultdict(set)
    extension_regexp = {".gs", ".vs", ".fs", ".shared", ".tsc", ".tse"}

    # get all *.gs *.vs *.fs *.shared from the two input directories
    shader_files = set(
        path
        for path in chain(shader_dir.glob("**/*"), shader_dir2.glob("**/*"))
        if path.suffix in extension_regexp
    )

    with (
        open(output_header, "w") as output_header_file,
        open(output_source, "w") as output_source_file,
    ):
        output_header_file.write(f"extern const char * {varname}[];\n")
        output_source_file.write(f"const char * {varname}[] = {{\n")

        for shader_file in shader_files:
            output_source_file.write(f'"{shader_file.name}", ""\n')

            contents = shader_file.read_text()
            for line in contents.splitlines():
                line = line.strip()

                # skip blank lines and obvious comments
                if not line or line.startswith("//") and not "*/" in line:
                    continue

                # write line, quoted, escaped and with a line feed
                escaped_line = line.replace("\\", "\\\\").replace('"', r"\"")
                output_source_file.write(f'"{escaped_line}\\n"\n')

                # include and ifdef dependencies
                if line.startswith("#include"):
                    include_deps[line.split()[1]].add(shader_file.name)
                elif line.startswith("#ifdef") or line.startswith("#ifndef"):
                    ifdef_deps[line.split()[1]].add(shader_file.name)

            output_source_file.write(",\n")
        output_source_file.write("0};\n")

        # include and ifdef dependencies
        for varname, deps in [
            ("_include_deps", include_deps),
            ("_ifdef_deps", ifdef_deps),
        ]:
            output_header_file.write(f"extern const char * {varname}[];\n")
            output_source_file.write(f"const char * {varname}[] = {{\n")
            for name, item_deps in deps.items():
                item_deps = '", "'.join(sorted(item_deps))
                output_source_file.write(f'"{name}", "{item_deps}", 0,\n')
            output_source_file.write("0};\n")


def create_buildinfo(output_dir_path: str, pymoldir: str = ".") -> None:
    output_dir = Path(output_dir_path)
    sha_raw = Popen(["git", "rev-parse", "HEAD"], cwd=pymoldir, stdout=PIPE).stdout
    sha = sha_raw.read().strip().decode() if sha_raw is not None else ""

    info_file = output_dir / "PyMOLBuildInfo.h"
    info_file.write_text(
        f"""
    #define _PyMOL_BUILD_DATE {time.time()}
    #define _PYMOL_BUILD_GIT_SHA "{sha}"
    """
    )


# handle extra arguments
def str2bool(v: str) -> bool:
    if v.lower() in ("true", "yes"):
        return True
    elif v.lower() in ("false", "no"):
        return False
    else:
        raise argparse.ArgumentTypeError("Boolean value expected.")


class options:
    osx_frameworks = True
    jobs = int(os.getenv("JOBS", 0))
    libxml = True
    glut = False
    use_msgpackc = "guess"
    testing = False
    openvr = False
    use_openmp = "no" if MAC else "yes"
    use_vtkm = "no"
    vmd_plugins = True


parser = argparse.ArgumentParser()
parser.add_argument(
    "--glut", dest="glut", type=str2bool, help="link with GLUT (legacy GUI)"
)
parser.add_argument(
    "--osx-frameworks",
    dest="osx_frameworks",
    help="on MacOS use XQuartz instead of native frameworks",
    type=str2bool,
)
parser.add_argument(
    "--jobs",
    "-j",
    type=int,
    help="for parallel builds (defaults to number of processors)",
)
parser.add_argument(
    "--libxml",
    type=str2bool,
    help="skip libxml2 dependency, disables COLLADA export",
)
parser.add_argument("--use-openmp", type=str2bool, help="Use OpenMP")
parser.add_argument(
    "--use-vtkm",
    choices=("2.0", "2.1", "2.2", "2.3", "no"),
    help="Use VTK-m for isosurface generation",
)
parser.add_argument(
    "--use-msgpackc",
    choices=("c++11", "c", "guess", "no"),
    help="c++11: use msgpack-c header-only library; c: link against "
    "shared library; no: disable fast MMTF load support",
)
parser.add_argument("--testing", type=str2bool, help="Build C-level tests")
parser.add_argument("--openvr", dest="openvr", type=str2bool)
parser.add_argument(
    "--vmd-plugins",
    dest="vmd_plugins",
    type=str2bool,
    help="Disable VMD molfile plugins (libnetcdf dependency)",
)
options, sys.argv[1:] = parser.parse_known_args(namespace=options)


def get_prefix_path() -> list[str]:
    """
    Return a list of paths which will be searched for "include",
    "include/freetype2", "lib", "lib64" etc.
    """
    paths = []

    if (prefix_path := os.environ.get("PREFIX_PATH")) is not None:
        paths += prefix_path.split(os.pathsep)

    if sys.platform.startswith("freebsd"):
        paths += ["/usr/local"]

    if not options.osx_frameworks:
        paths += ["usr/X11"]

    if MAC:
        for prefix in ["/sw", "/opt/local", "/usr/local"]:
            if sys.base_prefix.startswith(prefix):
                paths += [prefix]

    if is_conda_env():
        if WIN:
            if "CONDA_PREFIX" in os.environ:
                paths += [os.path.join(os.environ["CONDA_PREFIX"], "Library")]
            paths += [os.path.join(sys.prefix, "Library")]

        paths += [sys.prefix] + paths

    paths += ["/usr"]

    return paths


def is_conda_env():
    return (
        "conda" in sys.prefix
        or "conda" in sys.version
        or "Continuum" in sys.version
        or sys.prefix == os.getenv("CONDA_PREFIX")
    )


def guess_msgpackc():
    for prefix in prefix_path:
        for suffix in ["h", "hpp"]:
            f = os.path.join(prefix, "include", "msgpack", f"version_master.{suffix}")

            try:
                m = re.search(r"MSGPACK_VERSION_MAJOR\s+(\d+)", open(f).read())
            except EnvironmentError:
                continue

            if m is not None:
                major = int(m.group(1))
                if major > 1:
                    return "c++11"

    return "no"


class CMakeExtension(Extension):
    def __init__(
        self,
        name,
        sources,
        include_dirs=[],
        libraries=[],
        library_dirs=[],
        define_macros=[],
        extra_link_args=[],
        extra_compile_args=[],
    ):
        # don't invoke the original build_ext for this special extension
        super().__init__(name, sources=[])
        self.sources = sources
        self.include_dirs = include_dirs
        self.libraries = libraries
        self.library_dirs = library_dirs
        self.define_macros = define_macros
        self.extra_link_args = extra_link_args
        self.extra_compile_args = extra_compile_args


class build_ext_pymol(build_ext):
    def initialize_options(self) -> None:
        super().initialize_options()
        if DEBUG and not WIN:
            self.debug = False

    def run(self):
        for ext in self.extensions:
            self.build_cmake(ext)

    def build_cmake(self, ext):
        cwd = Path().absolute()

        # these dirs will be created in build_py, so if you don't have
        # any python sources to bundle, the dirs will be missing
        name_split = ext.name.split(".")
        target_name = name_split[-1]
        build_temp = Path(self.build_temp) / target_name
        build_temp.mkdir(parents=True, exist_ok=True)
        extdir = Path(self.get_ext_fullpath(ext.name))

        extdir.parent.mkdir(parents=True, exist_ok=True)

        def concat_paths(paths):
            return "".join(path.replace("\\", "/") + ";" for path in paths)

        config = "Debug" if DEBUG else "Release"
        all_files = ext.sources
        all_src = concat_paths(all_files)
        all_defs = "".join(mac[0] + ";" for mac in ext.define_macros)
        all_libs = "".join(f"{lib};" for lib in ext.libraries)
        all_ext_link = " ".join(ext.extra_link_args)
        all_comp_args = "".join(f"{arg};" for arg in ext.extra_compile_args)
        all_lib_dirs = concat_paths(ext.library_dirs)
        all_inc_dirs = concat_paths(ext.include_dirs)

        cmake_args = [
            f"-DTARGET_NAME={target_name}",
            f"-DCMAKE_BUILD_TYPE={config}",
            f"-DALL_INC_DIR={all_inc_dirs}",
            f"-DALL_SRC={all_src}",
            f"-DALL_DEF={all_defs}",
            f"-DALL_LIB_DIR={all_lib_dirs}",
            f"-DALL_LIB={all_libs}",
            f"-DALL_COMP_ARGS={all_comp_args}",
            f"-DALL_EXT_LINK={all_ext_link}",
        ]

        # example of build args
        build_args = ["--config", config]
        if not WIN:  # Win /MP flag on compilation level
            cpu_count = os.cpu_count() or 1
            build_args += [f"-j{cpu_count}"]

        os.chdir(str(build_temp))
        self.spawn(["cmake", str(cwd)] + cmake_args)
        if not self.dry_run:
            self.spawn(["cmake", "--build", "."] + build_args)

        os.chdir(str(cwd))


class build_py_pymol(build_py):
    def run(self):
        build_py.run(self)


class install_pymol(install):
    pymol_path = None
    bundled_pmw = False
    no_launcher = False

    user_options = install.user_options + [
        ("pymol-path=", None, "PYMOL_PATH"),
        ("bundled-pmw", None, "install bundled Pmw module"),
        ("no-launcher", None, "skip installation of the pymol launcher"),
    ]

    def finalize_options(self):
        install.finalize_options(self)

        self.pymol_path_is_default = self.pymol_path is None

        if self.pymol_path is None:
            self.pymol_path = os.path.join(self.install_libbase, "pymol", "pymol_path")
        elif self.root is not None:
            self.pymol_path = install_pymol.change_root(self.root, self.pymol_path)

    def run(self):
        super().run()

        assert self.pymol_path is not None
        self.install_pymol_path(self.pymol_path)

        if not self.no_launcher:
            self.make_launch_script()

        if self.bundled_pmw:
            raise Exception(
                "--bundled-pmw has been removed, please install Pmw from "
                "https://github.com/schrodinger/pmw-patched"
            )

    def unchroot(self, name):
        if self.root is not None and name.startswith(self.root):
            return name[len(self.root) :]
        return name

    def copy_tree_nosvn(self, src, dst):
        def ignore(src, names):
            return set([]).intersection(names)

        if os.path.exists(dst):
            shutil.rmtree(dst)
        print("copying %s -> %s" % (src, dst))
        shutil.copytree(src, dst, ignore=ignore)

    def copy(self, src, dst):
        copy = self.copy_tree_nosvn if os.path.isdir(src) else self.copy_file
        copy(src, dst)

    def install_pymol_path(self, base_path):
        self.mkpath(base_path)
        for name in [
            "LICENSE",
            "data",
            "test",
            "examples",
        ]:
            self.copy(name, os.path.join(base_path, name))

        if options.openvr:
            self.copy("contrib/vr/README.md", os.path.join(base_path, "README-VR.txt"))

    def make_launch_script(self):
        if sys.platform.startswith("win"):
            launch_script = "pymol.bat"
        else:
            launch_script = "pymol"

        self.mkpath(self.install_scripts)
        launch_script = os.path.join(self.install_scripts, launch_script)

        python_exe = os.path.abspath(sys.executable)
        site_packages_dir = sysconfig.get_path("purelib")
        pymol_file = self.unchroot(
            os.path.join(site_packages_dir, "pymol", "__init__.py")
        )
        pymol_path = self.unchroot(self.pymol_path)

        with open(launch_script, "w") as out:
            if WIN:
                # paths relative to launcher, if possible
                try:
                    python_exe = "%~dp0\\" + os.path.relpath(
                        python_exe, self.install_scripts
                    )
                except ValueError:
                    pass
                try:
                    pymol_file = "%~dp0\\" + os.path.relpath(
                        pymol_file, self.install_scripts
                    )
                except ValueError:
                    pymol_file = os.path.abspath(pymol_file)

                if not self.pymol_path_is_default:
                    out.write(f"set PYMOL_PATH={pymol_path}" + os.linesep)
                out.write('"%s" "%s"' % (python_exe, pymol_file))
                out.write(" %*" + os.linesep)
            else:
                out.write("#!/bin/sh" + os.linesep)
                if not self.pymol_path_is_default:
                    out.write(f'export PYMOL_PATH="{pymol_path}"' + os.linesep)
                out.write('exec "%s" "%s" "$@"' % (python_exe, pymol_file) + os.linesep)

        os.chmod(launch_script, 0o755)


# ============================================================================


# should be something like (build_base + "/generated"), but that's only
# known to build and install instances
generated_dir = os.path.join(os.environ.get("PYMOL_BLD", "build"), "generated")

create_all(generated_dir)

# can be changed with environment variable PREFIX_PATH
prefix_path = get_prefix_path()


pymol_src_dirs = [
    "ov/src",
    "layer0",
    "layer1",
    "layer2",
    "layer3",
    "layer4",
    "layer5",
    generated_dir,
]
libs = ["png", "freetype"]

inc_dirs = []
def_macros = []
lib_dirs = []
ext_comp_args = []
ext_link_args = []
ext_objects = []
data_files = []
ext_modules = []

if options.use_openmp == "yes":
    def_macros += [
        ("PYMOL_OPENMP", None),
    ]
    if MAC:
        ext_comp_args += ["-Xpreprocessor", "-fopenmp"]
        libs += ["omp"]
    elif WIN:
        ext_comp_args += ["/openmp"]
    else:
        ext_comp_args += ["-fopenmp"]
        ext_link_args += ["-fopenmp"]

if options.libxml:
    # COLLADA support
    def_macros += [("_HAVE_LIBXML", None)]
    libs += ["xml2"]

if options.use_msgpackc == "guess":
    options.use_msgpackc = guess_msgpackc()

if options.use_msgpackc == "no":
    def_macros += [("_PYMOL_NO_MSGPACKC", None)]
else:
    if options.use_msgpackc == "c++11":
        def_macros += [
            ("MMTF_MSGPACK_USE_CPP11", None),
            ("MSGPACK_NO_BOOST", None),
        ]
    else:
        libs += ["msgpackc"]

    pymol_src_dirs += ["contrib/mmtf-c"]

if not options.glut:
    def_macros += [
        ("_PYMOL_NO_MAIN", None),
    ]

if options.testing:
    pymol_src_dirs += ["layerCTest"]
    def_macros += [("_PYMOL_CTEST", None)]

if options.openvr:
    def_macros += [("_PYMOL_OPENVR", None)]
    pymol_src_dirs += [
        "contrib/vr",
    ]

inc_dirs += pymol_src_dirs

# ============================================================================
if MAC:
    libs += ["GLEW"]

    if options.osx_frameworks:
        ext_link_args += [
            "-framework OpenGL",
        ] + (options.glut) * [
            "-framework GLUT",
        ]
        def_macros += [
            ("_PYMOL_OSX", None),
        ]
    else:
        libs += [
            "GL",
        ] + (options.glut) * [
            "glut",
        ]

if WIN:
    libs.clear()

    libs += [
        "Advapi32",  # Registry (RegCloseKey etc.)
        "Ws2_32",  # htonl
    ]

    libs += (
        [
            "glew32",
            "freetype",
            "libpng",
        ]
        + (options.glut)
        * [
            "freeglut",
        ]
        + (options.libxml)
        * [
            "libxml2",
        ]
    )

    libs += [
        "opengl32",
    ]

if not (MAC or WIN):
    libs += [
        "GL",
        "GLEW",
    ] + (options.glut) * [
        "glut",
    ]

if options.use_vtkm != "no":
    for prefix in prefix_path:
        vtkm_inc_dir = os.path.join(prefix, "include", f"vtkm-{options.use_vtkm}")
        if os.path.exists(vtkm_inc_dir):
            break
    else:
        raise LookupError(
            f"VTK-m headers not found. PREFIX_PATH={':'.join(prefix_path)}"
        )
    def_macros += [
        ("_PYMOL_VTKM", None),
    ]
    inc_dirs += [
        vtkm_inc_dir,
        vtkm_inc_dir + "/vtkm/thirdparty/diy/vtkmdiy/include",
        vtkm_inc_dir + "/vtkm/thirdparty/lcl/vtkmlcl",
    ]
    libs += [
        f"vtkm_cont-{options.use_vtkm}",
        f"vtkm_filter_contour-{options.use_vtkm}",
        f"vtkm_filter_core-{options.use_vtkm}",
    ]

if options.vmd_plugins:
    # VMD plugin support
    inc_dirs += [
        "contrib/uiuc/plugins/include",
    ]
    pymol_src_dirs += [
        "contrib/uiuc/plugins/molfile_plugin/src",
    ]
    def_macros += [
        ("_PYMOL_VMD_PLUGINS", None),
    ]
    libs += [
        "netcdf",
    ]

if options.openvr:
    libs += [
        "openvr_api",
    ]

for prefix in prefix_path:
    for dirs, suffixes in [
        [
            inc_dirs,
            [
                ("include",),
                ("include", "freetype2"),
                ("include", "libxml2"),
                ("include", "openvr"),
            ],
        ],
        [lib_dirs, [("lib64",), ("lib",)]],
    ]:
        dirs.extend(filter(os.path.isdir, [os.path.join(prefix, *s) for s in suffixes]))


def get_pymol_version() -> str:
    return re.findall(r'_PyMOL_VERSION "(.*)"', open("layer0/Version.h").read())[0]


def get_sources(subdirs, suffixes=(".c", ".cpp")) -> list[str]:
    return sorted(
        [f for d in subdirs for s in suffixes for f in glob.glob(d + "/*" + s)]
    )


ext_modules += [
    CMakeExtension(
        name="_cmd",
        sources=get_sources(pymol_src_dirs),
        include_dirs=inc_dirs,
        libraries=libs,
        library_dirs=lib_dirs,
        define_macros=def_macros,
        extra_link_args=ext_link_args,
        extra_compile_args=ext_comp_args,
    ),
]

setup(
    cmdclass={
        "build_ext": build_ext_pymol,
        "build_py": build_py_pymol,
        "install": install_pymol,
    },
    version=get_pymol_version(),
    ext_modules=ext_modules,
)
