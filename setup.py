from distutils.file_util import copy_file
import os, glob

__author__ = 'marius'

pjoin = os.path.join

snufflings = pjoin(os.getenv('HOME'), '.snufflings')
cwd = os.getcwd()
files = glob.glob(pjoin(cwd, 'derec/*'))

if not os.path.exists(pjoin(snufflings, 'derec')):
    os.makedirs(pjoin(snufflings, 'derec'))

for fn in files:
    copy_file(src=os.path.join(cwd, 'derec', fn.rsplit('/').pop()),
              dst=os.path.join(snufflings, 'derec',fn.rsplit('/').pop()),
              link='sym')
