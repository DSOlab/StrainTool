import setuptools

setuptools.setup(name='pystrain',
      version='1.0-r1',
      description='Python Strain Tensor estimation tool.',
      url='https://github.com/DSOlab/StrainTool.git',
      author='Xanthos Papanikolaou, Dimitris Anastasiou',
      author_email='xanthos@mail.ntua.gr, dganastasiou@gmail.com',
      packages=setuptools.find_packages(),#['pystrain', 'pystrain.geodesy', 'pystrain.iotools'],
      install_requires=['numpy', 'scipy']
      )
