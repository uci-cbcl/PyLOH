from distutils.core import setup

setup(
      name = 'PyLOH',
      version = '1.0',
      description = 'Estimating tumor cellular frequency through copy number variation and loss of heterozygosity of cancer sequencing data.',
      author = 'Yi Li',
      author_email = 'yil8@uci.edu',
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