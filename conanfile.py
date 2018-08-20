import os
from conans import ConanFile, CMake
from conans.errors import ConanException


class XmsinterpConan(ConanFile):
    name = "xmsmesh"
    version = None
    license = "XMSNG Software License"
    url = "https://github.com/Aquaveo/xmsmesh"
    description = "Meshing library for XMS products"
    settings = "os", "compiler", "build_type", "arch"
    options = {
        "xms": [True, False],
        "pybind": [True, False],
        "testing": [True, False],
    }
    default_options = "xms=False", "pybind=False", "testing=False"
    generators = "cmake"
    build_requires = "cxxtest/4.4@aquaveo/stable"
    exports = "CMakeLists.txt", "LICENSE", "test_files/*"
    exports_sources = "xmsmesh/*"

    def configure(self):
        # Set verion dynamically using XMS_VERSION env variable.
        self.version = self.env.get('XMS_VERSION', 'master')

        # Raise ConanExceptions for Unsupported Versions
        s_os = self.settings.os
        s_compiler = self.settings.compiler
        s_compiler_version = self.settings.compiler.version

        self.options['xmscore'].xms = self.options.xms
        self.options['xmscore'].pybind = self.options.pybind
        self.options['xmscore'].testing = self.options.testing

        self.options['xmsinterp'].xms = self.options.xms
        self.options['xmsinterp'].pybind = self.options.pybind
        self.options['xmsinterp'].testing = self.options.testing

        if s_compiler == "apple-clang" and s_os == 'Linux':
            raise ConanException("Clang on Linux is not supported.")

        if s_compiler == "apple-clang" \
                and s_os == 'Macos' \
                and float(s_compiler_version.value) < 9.0:
            raise ConanException("Clang > 9.0 is required for Mac.")

    def requirements(self):
        """Requirments"""
        # If building for XMS, use the older, custom boost
        if self.options.xms:
            self.requires("boost/1.60.0@aquaveo/testing")
        else:
            self.requires("boost/1.66.0@conan/stable")

        # Pybind if not visual studio 2013
        if not (self.settings.compiler == 'Visual Studio' \
                and self.settings.compiler.version == "12") \
                and self.options.pybind:
            self.requires("pybind11/2.2.2@aquaveo/stable")

        self.requires("xmscore/[>=1.0.33]@aquaveo/stable")
        self.requires("xmsinterp/[>1.0.7]@aquaveo/stable")

    def build(self):
        cmake = CMake(self)

        if self.settings.compiler == 'Visual Studio' \
           and self.settings.compiler.version == "12":
            cmake.definitions["XMS_BUILD"] = self.options.xms

        # CXXTest doesn't play nice with PyBind. Also, it would be nice to not
        # have tests in release code. Thus, if we want to run tests, we will
        # build a test version (without python), run the tests, and then (on
        # sucess) rebuild the library without tests.
        cmake.definitions["IS_PYTHON_BUILD"] = self.options.pybind
        cmake.definitions["BUILD_TESTING"] = self.options.testing
        cmake.definitions["PYTHON_TARGET_VERSION"] = self.env.get("PYTHON_TARGET_VERSION", "3.6")
        cmake.configure(source_folder=".")
        cmake.build()

        if self.options.testing:
            print("***********(0.0)*************")
            try:
                cmake.test()
            except ConanException:
                raise
            finally:
                if os.path.isfile("TEST-cxxtest.xml"):
                    with open("TEST-cxxtest.xml", "r") as f:
                        for line in f.readlines():
                            no_newline = line.strip('\n')
                            print(no_newline)
                print("***********(0.0)*************")

    def package(self):
        self.copy("*.h", dst="include/xmsmesh", src="xmsmesh")
        self.copy("*.lib", dst="lib", keep_path=False)
        self.copy("*.dll", dst="bin", keep_path=False)
        self.copy("*.dylib*", dst="lib", keep_path=False)
        self.copy("*.so", dst="lib", keep_path=False)
        self.copy("*.a", dst="lib", keep_path=False)
        self.copy("license", dst="licenses", ignore_case=True, keep_path=False)

    def package_info(self):
        if self.settings.build_type == 'Debug':
            self.cpp_info.libs = ["xmsmesh_d"]
        else:
            self.cpp_info.libs = ["xmsmesh"]
