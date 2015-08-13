import os
import paths
import glob

figstr = """
\Figure{{figures/pointsource_seds/{fname}}}
{{Point source photometry of {objname}.  See Figure \\ref{{fig:d4sed}} for
details}}
{{fig:{objname}sed}}{{1}}{{7in}}\n\clearpage\n"""

with open(paths.texpath('ptsrcphotometryfigures.tex'),'w') as f:
    for fn in glob.glob(paths.fpath('pointsource_seds/*png')):
        fname = os.path.split(fn)[-1]
        objname = fname.split("_")[0]
        f.write(figstr.format(fname=fname, objname=objname))


