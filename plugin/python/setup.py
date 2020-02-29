from distutils.core import setup
from distutils.extension import Extension
import os
import sys
import platform

openmm_dir = '@OPENMM_DIR@'
seekrplugin_header_dir = '@SEEKRPLUGIN_HEADER_DIR@'
seekrplugin_library_dir = '@SEEKRPLUGIN_LIBRARY_DIR@'

# setup extra compile and link arguments on Mac
extra_compile_args = ['-std=c++11']
extra_link_args = []

if platform.system() == 'Darwin':
    extra_compile_args += ['-stdlib=libc++', '-mmacosx-version-min=10.7']
    extra_link_args += ['-stdlib=libc++', '-mmacosx-version-min=10.7', '-Wl', '-rpath', openmm_dir+'/lib']

extension = Extension(name='_seekrplugin',
                      sources=['SeekrPluginWrapper.cpp'],
                      libraries=['OpenMM', 'SeekrPlugin'],
                      include_dirs=[os.path.join(openmm_dir, 'include'), seekrplugin_header_dir],
                      library_dirs=[os.path.join(openmm_dir, 'lib'), seekrplugin_library_dir],
                      extra_compile_args=extra_compile_args,
                      extra_link_args=extra_link_args
                     )

setup(name='seekrplugin',
      version='0.1',
      py_modules=['seekrplugin'],
      ext_modules=[extension],
     )
