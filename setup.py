from distutils.core import setup

setup(
      name = 'PyLOH',
      version = '1.1',
      description = 'Deconvolving tumor purity and ploidy by integrating copy number alterations and loss of heterozygosity',
      author = 'Yi Li, Andrew Roth',
      author_email = 'yil8@uci.edu, andrewjlroth@gmail.com',
      url = 'https://github.com/uci-cbcl/PyLOH',
      license = 'GNU GPL v2',
      packages = [
                  'pyloh',
                  'pyloh.preprocess',
                  'pyloh.model',
                  'pyloh.postprocess'
      ],
      scripts = ['PyLOH.py']
)
