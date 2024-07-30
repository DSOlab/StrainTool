import setuptools

setuptools.setup(name='pystrain',
      version='1.1beta1',
      description='Python Strain Tensor estimation tool.',
      url='https://github.com/DSOlab/StrainTool.git',
      author='Xanthos Papanikolaou, Dimitris Anastasiou',
      author_email='xanthos@mail.ntua.gr, danastasiou@mail.ntua.gr',
      packages=setuptools.find_packages(),#['pystrain', 'pystrain.geodesy', 'pystrain.iotools'],
      install_requires=['numpy', 'scipy', 'utm']
      )
