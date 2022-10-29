from setuptools import setup

def readme():
    with open('README.md') as f:
        return f.read()

setup(name='gssnng',
      version='0.2.1',
      description='Gene Set Scoring on the Nearest Neighbor Graph (gssnng)',
      url='http://github.com/gibbsdavidl/gssnng',
      author='David Gibbs,Michael Strasser',
      author_email='gibbsdavidl@gmail.com',
      license='MIT',
      packages=['gssnng'],
      install_requires=[
          'pandas', 'numpy', 'matplotlib', 'seaborn', 'scipy', 'statsmodels', 'scanpy', 'tqdm'
      ],
      zip_safe=False)
