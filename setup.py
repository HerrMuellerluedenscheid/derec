from distutils.file_util import copy_file
from distutils.core import setup, Command
import os, glob

__author__ = 'marius'

pjoin = os.path.join


class SetupBuildCommand(Command):
    """
    Master setup build command to subclass from.
    """

    user_options = []

    def initialize_options(self):
        """
        Setup the current dir.
        """
        self._dir = os.getcwd()

    def finalize_options(self):
        """
        Set final values for all the options that this command supports.
        """
        pass


class LinkSnufflingFiles(SetupBuildCommand):

    description = "Create symbolic links and subdirectory in $HOME/.snufflings"

    def run(self):
        snufflings = pjoin(os.getenv('HOME'), '.snufflings')
        cwd = os.getcwd()
        files = glob.glob(pjoin(cwd, 'derec/*'))

        if not os.path.exists(pjoin(snufflings, 'derec')):
            os.makedirs(pjoin(snufflings, 'derec'))

        for fn in files:
            copy_file(src=os.path.join(cwd, 'derec', fn.rsplit('/').pop()),
                      dst=os.path.join(snufflings, 'derec',fn.rsplit('/').pop()),
                      link='sym')


setup(name='derec',
      version='1.0',
      description='depth relocation',
      packages=['derec'],

      cmdclass={'link': LinkSnufflingFiles},
      )
