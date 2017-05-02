from setuptools import setup

setup(name='pystrain',
      version='0.1',
      description='Python Strain Tensor estimation tool.',
      url='https://github.com/DSOlab/StrainTool.git',
      author='Xanthos, Mitsos',
      author_email='xanthos@mail.ntua.gr, danast@mail.ntua.gr',
      license='MIT',
      packages=['pystrain', 'pystrain.geodesy', 'pystrain.iotools'],
      install_requires=[],
      zip_safe=False)
