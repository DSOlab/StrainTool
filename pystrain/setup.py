import setuptools

setuptools.setup(name='pystrain',
      version='1.0',
      description='Python Strain Tensor estimation tool.',
      url='https://github.com/DSOlab/StrainTool.git',
      author='Xanthos Papanikolaou, Demitris Anastasiou',
      author_email='xanthos@mail.ntua.gr, danast@mail.ntua.gr',
      packages=setuptools.find_packages(),#['pystrain', 'pystrain.geodesy', 'pystrain.iotools'],
      install_requires=['numpy', 'scipy']
      )
