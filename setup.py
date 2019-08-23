from distutils.core import setup, Extension

setup(name='airmass',
      install_requires=['numpy', 'pandas'],
      version='1.0',
      description='A Python module for computing airmass on Earth, originally from Reed D. Meyer.',
      ext_modules=[Extension('airmassc', sources = ['extensions/airmassc.c'])])
