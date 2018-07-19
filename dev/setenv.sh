#!/usr/bin/env bash
export XMS_VERSION="0.0.0"
export CONAN_REFERENCE="xmsmesh/${XMS_VERSION}"
export CONAN_USERNAME="aquaveo"
export CONAN_CHANNEL="stable"
export CONAN_GCC_VERSIONS="5"
export CONAN_ARCHS="x86_64"
export CONAN_BUILD_TYPES="Debug"
#export CONAN_UPLOAD="1"
#export CONAN_DOCKER_IMAGE="lasote/conangcc${CONAN_GCC_VERSIONS}"
printenv | grep XMS
printenv | grep CONAN
