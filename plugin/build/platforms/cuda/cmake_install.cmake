# Install script for directory: /home/lvotapka/openseekr/plugin/platforms/cuda

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/home/lvotapka/bin/openmm")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  if(EXISTS "$ENV{DESTDIR}/home/lvotapka/bin/openmm/lib/plugins/libSeekrPluginCUDA.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/home/lvotapka/bin/openmm/lib/plugins/libSeekrPluginCUDA.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}/home/lvotapka/bin/openmm/lib/plugins/libSeekrPluginCUDA.so"
         RPATH "")
  endif()
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/lvotapka/bin/openmm/lib/plugins/libSeekrPluginCUDA.so")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/home/lvotapka/bin/openmm/lib/plugins" TYPE SHARED_LIBRARY FILES "/home/lvotapka/openseekr/plugin/build/platforms/cuda/libSeekrPluginCUDA.so")
  if(EXISTS "$ENV{DESTDIR}/home/lvotapka/bin/openmm/lib/plugins/libSeekrPluginCUDA.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/home/lvotapka/bin/openmm/lib/plugins/libSeekrPluginCUDA.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}/home/lvotapka/bin/openmm/lib/plugins/libSeekrPluginCUDA.so"
         OLD_RPATH "/home/lvotapka/openseekr/plugin/build:/home/lvotapka/bin/openmm/lib:/home/lvotapka/bin/openmm/lib/plugins:"
         NEW_RPATH "")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}/home/lvotapka/bin/openmm/lib/plugins/libSeekrPluginCUDA.so")
    endif()
  endif()
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/home/lvotapka/openseekr/plugin/build/platforms/cuda/tests/cmake_install.cmake")

endif()

