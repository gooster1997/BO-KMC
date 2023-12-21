from distutils.core import setup, Extension

setup(name='KMC',
      ext_modules=[
        Extension('KMC',
                  ['monte_carlo.c'],
                  )
        ]
)