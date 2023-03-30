from distutils.core import setup, Extension

def main():
    setup(name="Monte",
          version="1.0.0",
          description="Python interface for the Monte under C",
          author="gooster",
          ext_modules=[Extension("Monte", ["therm_coup.c", "therm_disp.c", "photo_coup.c", "photo_disp.c"])])
          
if __name__=="__main__":
    main()
          
