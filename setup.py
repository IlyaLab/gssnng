from setuptools import setup

def readme():
    with open('README.rst') as f:
        return f.read()

setup(name='nnggss',
      version='0.0.1',
      description='Nearest Neighbor Graph Gene Set Scoring (nnggss)',
      url='http://github.com/gibbsdavidl/nnggss',
      author='David Gibbs,Michael Strasser',
      author_email='gibbsdavidl@gmail.com',
      license='MIT',
      packages=['nnggss'],
      install_requires=[
          'pandas', 'numpy', 'matplotlib', 'seaborn', 'scipy', 'statsmodels', 'scanpy', 'tqdm'
      ],
      zip_safe=False)
