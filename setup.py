from distutils.core import setup

setup(name='checkrna',
      version='0.1.0',
      description='A tool to re-reference 13C chemical shifts of RNA',
      author='Alejandro Icazatti and Osvaldo Martin',
      author_email='ale.icazatti@gmail.com, aloctavodia@gmail.com',
      url='https://github.com/BIOS-IMASL/13Check_RNA',
      install_requires=['numpy','pandas','pynmrstar'],
      packages=['checkrna'],
     )
