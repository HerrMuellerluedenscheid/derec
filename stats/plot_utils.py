import re

def to_percent(y, position):
    # Ignore the passed in position. This has the effect of scaling the default
    # tick locations.
    s = str(100 * y)

    # The percent symbol needs escaping in latex
    if matplotlib.rcParams['text.usetex'] == True:
        return s + r'$\%$'
    else:
        return s + '%'

def isplain_dir(dirstr):
    return re.match('.*plain.*', dirstr)!=None
def iscastor_dir(dirstr):
    return re.match('.*castor.*', dirstr)!=None

def isdoctar_dir(dirstr):
    return re.match('.*doctar.*', dirstr)!=None

