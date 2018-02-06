#!/usr/bin/python

###  Copyright 2014-2018 Manuel N. Melo ###
#########################################################################
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.
#########################################################################


import sys
import argparse
import os
import MDAnalysis
import math
import numpy
import re
import datetime



# Static descriptions ##################################################
########################################################################

desc = '''
Counts crossing events across a finite width xy-plane membrane.
The system is assumed periodic in z but the pbc in the trajectory must have been
removed (trjconv option -pbc nojump).
Two index groups define the top and bottom of the membrane. A crossing is only
counted when an atom exits on the opposite side of the membrane to the one it
came in.
When analyzing a single channel, the option -mult can be passed to restrict
counting to the cylinder of radius mult*RoG, where RoG is the radius of gyration
(not mass weighted) of the delimiting groups. Care must be taken to ensure these
groups remain whole for the whole trajectory.
fluxer.py assumes an xy flat membrane slab. Other topologies or heavy
fluctuations will produce erroneous results.

The schematic of the counting is as follows:
   _ _ _ _ _ _ _ _ _____
  |##Slab 4####| |#|  B
  | ~ ~ ~ ~ ~ ~ ~ ~|  O
  |~ Slab 3  ~ ~ ~ |  X  
  | ~ ~ ~ ~ ~ ~ ~ ~|  1 
  |##Slab 2####| |#|___ 
  |##Slab 2####| |#|  B  LEGEND:
  | ~ ~ ~ ~ ~ ~ ~ ~|  O   ~ = Water
  |~ Slab 1  ~ ~E~ |  X   # = Membrane                            
  | ~ ~ ~ ~ N ~ ~ ~|  0  | |= Channel cylinder                   
  |##Slab 0##n#| |#|___   E = Eligible molecule, in water        
  |##Slab 0####|e|#|  B   N = Not eligible molecule, in water    
  | ~ ~ ~ ~ ~ ~ ~ ~|  O   e = Eligible molecule, in membrane     
  |~ Slab -1 ~ ~ ~ |  X   n = not eligible molecule, in membrane 
  | ~ ~ ~ ~ ~ ~ ~ ~| -1 
  |##Slab -2###| |#|___ 
  |##Slab -2###| |#|  B 
  | ~ ~ ~ ~ ~ ~ ~ ~|  O 
  |~ Slab -3 ~ ~ ~ |  X 
  | ~ ~ ~ ~ ~ ~ ~ ~| -2
  |##Slab -4###| |#|___ 

An event is recorded if a water molecule arrives at the water slab of a box,
having last been at the water slab of an adjacent box. In addition, the crossing
must be in line with the channel, meaning the molecule must have been in the
channel ('eligible') in the frame preceding the one when it reached the
water slab.

Columns in the output .xvg and the analysis total are as follows:
1:Time
2:+Flux
3:-Flux
4:Crossings that were bigjumps
5:Total bigjumps
6-22:Jumptype discrimination.

Bigjump signals crossings between non-consecutive regions (namely, directly
jumping from one water slab to the next in one step). Having a large number of
these may indicate too large a frame step for the analysis. 300 ps per frame
is a typically safe value (for water).

The jumptypes code for different events that are recorded. These are numbered
0-16 and printed out in order in the .xvg columns 6-22.
Each frame the slab and eligibility of each molecule is compared to the previous
frame's. The jumptypes codes are as per the following matrix (the 'e E n N'
labels are as per the legend above):

         To(new):   
From:     e   E   n   N    +: stands for a clean crossing;                
(old) e *16 +15 *14 +13    *: is a crossing that can only happen
      E ?12 *11 ?10  *9       together with a bigjump;
      n  x8  x7  x6  x5    x: will not be counted as a crossing.
      N  ?4  *3  x2  x1    ?: is a crossing that should never be observed
                              (out of internal consistency - it requires
                              a jump so big that it is prevented by the
                              pbc treatment of nojump).

In addition to these, jumptype = 0 means no slab difference, or a membrane-water
transition back to the last water slab the molecule visited.
Eligibility is defined as being within the z-axis-oriented cylinder defined by
the ref groups' CoG and RoG (scaled by the -mult option). 
The assignment of '*' and 'x' is somewhat arbitrary regarding some exotic jumps
that will only be of concern when the frame step is so high that the results
will be rubbish anyway. The bottom line is that one should aim for mostly
'+' events; too many bigjumps indicate that the step is inappropriate.
   
'''

vernum = "v2013.08.15"

ver = '''%s
by Manuel Melo (m.n.melo@rug.nl)''' % (vernum)

header = '''#   Produced by fluxer.py %s (by Manuel Melo; m.n.melo@rug.nl)
#   %s
#\n''' % (vernum, " ".join(sys.argv))

fluxheader   = "#\n#   Flux vs time.\n"
fluxprologue = "# Columns: 1:Time, 2:Upward flux, 3:Downward flux, 4:Bigjumps counted as flux, 5:Total bigjumps, 6-22:Jumptype discrimination.\n"

cumfluxheader   = "#\n#   Cumulative flux vs time.\n"
cumfluxprologue = "# Columns: 1:Time, 2:Total cum. flux, 3:Upward cum. flux, 4:Downward cum. flux.\n"


# Main #################################################################
########################################################################

# flx_syst is a global instance of FluxSystem with all everyone needs to know about everyone.
# opts is the argparse options object.
flx_syst = None
opts = None
def main(argv=None):
    global opts
    global flx_syst

    ## Set up option checking
    def check_file(fname):
        if not os.path.exists(fname):
            sys.exit('Error: Can\'t find file %s' % (fname))
        if not os.access(fname, os.R_OK):
            sys.exit('Error: Permission denied to read file %s' % (fname))
        return fname
    def check_outfile(fname):
        dirname = os.path.dirname(fname)
        if not dirname:
            dirname = '.'
        if not os.access(dirname, os.W_OK):
            sys.exit('Error: Permission denied to write file %s' % (fname))
        return fname
    def check_positive(val):
        if val < 0:
            sys.exit('Error: Argument must be >= 0: %r' % (val))

    ## Options
    parser = argparse.ArgumentParser(description=desc, formatter_class=argparse.RawDescriptionHelpFormatter, epilog=" ", prog=argv[0])
    # Argparse is stupid and allows EITHER a properly formatted help section,
    #  or the addition of the defaults to the --help argument list. Not both...
    #parser = argparse.ArgumentParser(description=desc, formatter_class=argparse.ArgumentDefaultsHelpFormatter, epilog=" ", prog=argv[0])
    parser.add_argument('-f', metavar='TRAJ', dest='xtc', default='traj.xtc',
            help = 'file\tThe trajectory to analyze.')
    parser.add_argument('-c', metavar='STRUCT', dest='gro', default='confout.gro',
            help = 'file\t.gro or .pdb file with the same atom numbering as the trajectory.')
    parser.add_argument('-o', metavar='OUT', dest='outfile', default='flux.xvg',
            help = 'file\tThe flux vs time file. Read the --help for the meaning of each column in the output.')
    parser.add_argument('-oc', metavar='OUT', dest='cumfile', default='flux-cumul.xvg',
            help = 'file\tThe cumulative flux vs time file. Columns are: time, total flux, upward flux, downward flux.')
    parser.add_argument('-dt', metavar='STEP', type=float, dest='dt',
            help = 'real\tThe frame time step (DEPRECATED: default is now to do the right thing and read dt from trajectory).')
    parser.add_argument('-b', metavar='TIME', type=float, dest='starttime', default=0,
            help = 'real\tTime to begin analysis from.')
    parser.add_argument('-e', metavar='TIME', type=float, dest='endtime', default='inf',
            help = 'real\tTime to end analysis at.')
    parser.add_argument('-skip', metavar='FRAMES', type=int, dest='skip', default=1,
            help = 'int\tNumber of frames to skip when analyzing.')
    parser.add_argument('-atname', metavar='NAME', dest='watername', default='W',
            help = 'str\tThe atom name counted for the crossings.')
    parser.add_argument('-top', metavar='GRP', type=str.lower, dest='topgrp', default='top',
            help = 'str\tThe atom group name defining the membrane top.')
    parser.add_argument('-bottom', metavar='GRP', type=str.lower, dest='botgrp', default='bottom',
            help = 'str\tThe atom group name defining the membrane bottom.')
    parser.add_argument('-n', metavar='INDEX', dest='ndx', default='index.ndx',
            help = 'file\tIndex file where to read the membrane delimiting groups from.')
    parser.add_argument('-toroidal', metavar='THICKNESS', type=float, dest='toroidal', default = 0.,
            help = 'real\tAssumes a toroidal pore and sets its z-thickness, ignoring the other top/bottom limits. Only the \'top\' index group will be used and its CoM should then match the membrane\'s (use all the lipids, for example). Such channel-independent pores should always be analyzed with -mult set to 0, as no pore tracking is implemented (yet); no distinction will be made between through-the-membrane and through-the-pore crossings.')
    parser.add_argument('-noscale', action='store_false', dest='tor_scale',
            help = 'bool\tBy default the thickness defined with -toroidal is scaled with the box\'s Z dimension. Set to keep constant.')
    parser.add_argument('-mult', metavar='MULT', type=lambda x:float(x)**2, dest='radius_multsq', default=0,
            help = 'real\tMultiplier of the xy radius of gyration of the top and bottom groups to consider as part of the channel. Set to 0 (default) to count all membrane crossings.')
    parser.add_argument('-v', metavar='LEVEL', type=int, choices=[0,1,2], dest='verbose', default=1,
            help = 'enum\tVerbosity level. 0:quiet, 1:progress+totals, 2:debug')
    parser.add_argument('-noqc', action='store_false', dest='qc',
            help = 'bool\tBy default quality control is done by outputting individual counts of crossing types. Set to disable.')
    parser.add_argument('-nohdr', action='store_false', dest='xvgheaders',
            help = 'bool\tVersion info, used command line parameters, and column identification are added (#commented) to the header of the output .xvg files. Set to disable header inclusion.')
    parser.add_argument('-V', '--version', action='version', version='%%(prog)s %s'%ver,
            help = 'Prints the script version and exits.')
    # Parse
    opts = parser.parse_args(argv[1:])
    ## Post option handling
    check_file(opts.gro)
    check_file(opts.xtc)
    check_file(opts.ndx)
    check_outfile(opts.outfile)
    check_outfile(opts.cumfile)
    map(check_positive,(opts.starttime,opts.endtime,opts.skip,opts.toroidal))
    if opts.dt:
        check_positive(opts.dt)
    if opts.endtime < opts.starttime:
        sys.exit('Error: Endframe lower than starttime.')
    if len(opts.watername) > 5 or re.search(' ', opts.watername):
        sys.exit('Error: Invalid atom name %s' % opts.watername)
    if opts.toroidal and opts.radius_multsq:
        sys.exit('''Error: You are analyzing a toroidal pore and have set the -mult option. This is currently not implemented because the pore-delimiting atoms may change along the trajectory.
Use -mult 0 (the default) when using -toroidal. The downside is that all crossings will be counted, even those outside the pore, which will hopefully be negligible.''')

    # some more checking is done along the way.
    
    ## Let rip
    flx_syst=FluxSystem()
    flx_syst.go()

# Classes ##############################################################
########################################################################

class FluxSystem:
    def __init__(self):
        self.syst = MDAnalysis.Universe(opts.gro, opts.xtc)
        self.realdt = self.syst.trajectory.dt
        if opts.dt:
            self.dt = opts.dt
        else:
            self.dt = self.realdt

        self.rings = readndx([opts.topgrp,opts.botgrp], opts.ndx)
        if not len(self.rings[opts.topgrp]):
            sys.exit('Error: Empty or not-found index top group \'%s\'.' % opts.topgrp)
        if not (len(self.rings[opts.botgrp]) or opts.toroidal):
            sys.exit('Error: Empty or not-found index bottom group \'%s\'.' % opts.botgrp)
        self.topring_atgr = self.syst.atoms[self.rings[opts.topgrp]]
        self.botring_atgr = self.syst.atoms[self.rings[opts.botgrp]]
        self.bothrings_atgr = self.topring_atgr + self.botring_atgr
        self.ringatms_inv = 1./len(self.bothrings_atgr)
        self.waters = None
        if opts.toroidal and opts.tor_scale:
            self.rel_membrane_top = opts.toroidal/(2*self.syst.dimensions[2])
            self.rel_membrane_bot = -self.rel_membrane_top 

    def go(self):
        if opts.verbose:
            skipping = True
            sys.stderr.write("Skipping...\n")
            loop_time = ThenNow()
            loop_time.fill(datetime.datetime.now())
        for self.snapshot in self.syst.trajectory:
            if skipping and opts.verbose:
                sys.stderr.write("Frame %d  t= %.1f ps      \r" % (self.snapshot.frame, self.snapshot.time))

            if self.snapshot.time >= opts.starttime and not round((self.snapshot.time - opts.starttime)/self.dt) % opts.skip:
                if skipping and opts.verbose:
                    skipping = False
                    sys.stderr.write("\nAnalyzing...\n")

                if self.snapshot.time > opts.endtime:
                    break

                if opts.verbose:
                    loop_time.update(datetime.datetime.now())
                    loop_dtime = loop_time.new - loop_time.old
                    sys.stderr.write("Frame %d  t= %.1f ps  (%.3f s/frame)     \r" % (self.snapshot.frame, self.snapshot.time, float(loop_dtime.seconds)+loop_dtime.microseconds/1000000.))

                self.x_dim, self.y_dim, self.z_dim = self.snapshot.dimensions[:3]
                if opts.toroidal:
                    self.channel_cog = self.topring_atgr.center_of_geometry()
                    if not opts.tor_scale:
                        self.rel_membrane_top = opts.toroidal/(2*self.z_dim)
                        self.rel_membrane_bot = -self.rel_membrane_top 
                else:
                    self.channel_cog = (self.topring_atgr.center_of_geometry()+self.botring_atgr.center_of_geometry())/2
                    # xy-plane radius of gyration; not mass weighted
                    self.cutoffsq = opts.radius_multsq * self.ringatms_inv * numpy.sum(numpy.power((self.bothrings_atgr.positions-self.channel_cog)*[1,1,0], 2))
                    self.rel_membrane_top = (self.topring_atgr.center_of_geometry()[2]-self.channel_cog[2])/self.z_dim
                    self.rel_membrane_bot = -self.rel_membrane_top 

                if not self.waters:
                    self.waters = WaterCollection()
                    continue

                self.waters.update_flux()
                self.waters.output_frame()

        if opts.verbose==1:
            sys.stderr.write("\n")
        self.waters.output_close()

class ThenNow:
    def __init__(self, oldval=None, newval=None):
        self.set(oldval, newval)
    def set(self, oldval, newval):
        self.old = oldval
        self.new = newval
    def fill(self, val):
        # Fill variant for the initial case where we have to assign both at initialization.
        self.set(val, val)
    def update(self, val, fill=False):
        if fill:
            self.fill(val)
        else:
            self.old = self.new
            self.new = val
#
##
###
class WaterCollection:
#    #                                     Counted  Total    
#    #                         +Flux -Flux Bigjumps Bigjumps  
#    frame_tally = numpy.array([  0,    0,     0,      0   ])

    def __init__(self):
        try:
            self.water_atms = flx_syst.syst.select_atoms("name %s"% opts.watername)
        except MDAnalysis.NoDataError:
            sys.exit('Error: No atoms found with name \'%s\'.' % opts.watername)
        nwaters = len(self.water_atms)

        self.XVG = open(opts.outfile, 'w')
        self.CXVG = open(opts.cumfile, 'w')
        if opts.xvgheaders:
            self.XVG.write(fluxheader)
            self.XVG.write(header)
            self.XVG.write(fluxprologue)
            self.CXVG.write(cumfluxheader)
            self.CXVG.write(header)
            self.CXVG.write(cumfluxprologue)

        self.total_tally = numpy.array([0, 0, 0, 0 ]) 
        self.total_jtype = numpy.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
        self.jumptype_allow = numpy.array([0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1])
        self.inmembrane=ThenNow()
        self.slab=ThenNow()
        self.waterslab=ThenNow()
        self.waterslab.fill(numpy.zeros(nwaters, dtype='int'))
        self.eligible=ThenNow(True, True)
        self.inmembrane=ThenNow()
        self.initialize = True
        self.update_flux()
        self.initialize = False

    def update_flux(self):
        # bring in the centered coordinates
        self.coords = self.water_atms.positions - flx_syst.channel_cog

        # set the slab and waterslab
        rx = self.coords[:,2]/flx_syst.z_dim
        mdf = numpy.modf(rx)
        self.inmembrane.update(numpy.abs(mdf[0]-numpy.round(mdf[0])) < flx_syst.rel_membrane_top, self.initialize)
        self.slab.update((numpy.round(rx)*2*self.inmembrane.new + (mdf[1]*2+numpy.sign(rx))*(~self.inmembrane.new)).astype(numpy.int), self.initialize)
        slabdiff = self.slab.new-self.slab.old
        bigjump = numpy.abs(slabdiff) > 1
        wslab_updt0 = ~bigjump * self.inmembrane.new
        wslab_updt1 = bigjump * self.inmembrane.new
        self.waterslab.update(self.waterslab.new*wslab_updt0 + self.slab.new*(~self.inmembrane.new) + (self.slab.new-numpy.sign(slabdiff))*wslab_updt1, self.initialize)
        wslabdiff = self.waterslab.new - self.waterslab.old

        # set eligibility
        if opts.radius_multsq:
            x = self.coords[:,0]-numpy.round(self.coords[:,0]/flx_syst.x_dim)*flx_syst.x_dim
            y = self.coords[:,1]-numpy.round(self.coords[:,1]/flx_syst.y_dim)*flx_syst.y_dim
            self.eligible.update((numpy.power(x,2)+numpy.power(y,2)) < flx_syst.cutoffsq, self.initialize)

        # jumptype and flux
        jumptype =  (self.waterslab.old!=0) * (wslabdiff!=0) * (numpy.left_shift(self.eligible.old,3) + numpy.left_shift(self.inmembrane.old,2) +
                    numpy.left_shift(self.eligible.new, 1) + self.inmembrane.new + 1)
        flux = numpy.sign(wslabdiff)*self.jumptype_allow[jumptype]

        # tally up
        self.frame_jtype = numpy.bincount(jumptype, minlength = 17)
        #values: ((+Flux, -Flux, Bigjumps counted as jumps, Total bigjumps), jumptype)
        #  jumptype reports jump type, as per the from:to matrix above.
        #  jumptype = 0 means no slab difference (not the same as flux = 0 as
        #  there may be a slab difference but no jump, as in a
        #  notelig+even -> notelig+noteven transition).
        self.frame_tally = numpy.array([numpy.sum(flux>0), numpy.sum(flux<0), numpy.sum(bigjump), numpy.sum(numpy.multiply(bigjump,flux!=0))])
        self.total_tally += self.frame_tally
        self.CXVG.write("%f\t%d\t%d\t%d\n" % (flx_syst.snapshot.time*flx_syst.dt/flx_syst.realdt, self.total_tally[0]+self.total_tally[1],self.total_tally[0],self.total_tally[1]))
        self.total_jtype += self.frame_jtype

        if opts.verbose > 1:
            jump         = (self.waterslab.old!=0)*(wslabdiff!=0)
            # Add here water indices to be debugged every frame.
            for mol in []:
                jump[mol]    = True
            if numpy.sum(jump):
                sys.stderr.write("Frame: %d   %.1f ps  rel.top: %f  zdim: %f\n" % (flx_syst.snapshot.frame, flx_syst.snapshot.time*flx_syst.dt/flx_syst.realdt, flx_syst.rel_membrane_top, flx_syst.z_dim))
            dbgindices   = numpy.where(jump)
            dbgcoords    = self.coords.compress(jump, axis=0)
            dbginmemb_o  = self.inmembrane.old.compress(jump)
            dbginmemb_n  = self.inmembrane.new.compress(jump)
            dbgslab_o    = self.slab.old.compress(jump)
            dbgslab_n    = self.slab.new.compress(jump)
            dbgslab_d    = slabdiff.compress(jump)
            dbgwslab_o   = self.waterslab.old.compress(jump)
            dbgwslab_n   = self.waterslab.new.compress(jump)
            dbgwslab_d   = wslabdiff.compress(jump)
            dbgbjump     = bigjump.compress(jump)
            dbgjtype     = jumptype.compress(jump)

            for tgt_mol in range(len(dbgcoords)):
                sys.stderr.write("""W#: %d  coords: %r  inmembrane.o: %r  inmembrane.n: %r
   slab.o: %d   slab.n: %d   slabdiff: %d
   wslab.o: %d  wslab.n: %d  wslabdiff: %d
   bigjump: %r  jumptype: %d\n""" % (dbgindices[0][tgt_mol], dbgcoords[tgt_mol], dbginmemb_o[tgt_mol], dbginmemb_n[tgt_mol], dbgslab_o[tgt_mol], dbgslab_n[tgt_mol], dbgslab_d[tgt_mol], dbgwslab_o[tgt_mol], dbgwslab_n[tgt_mol], dbgwslab_d[tgt_mol], dbgbjump[tgt_mol], dbgjtype[tgt_mol])
)

    def output_frame(self):
        self.XVG.write("%f\t" % (flx_syst.snapshot.time*flx_syst.dt/flx_syst.realdt))
        self.XVG.write(("%d\t"*4) % tuple(self.frame_tally))
        if opts.qc:
            self.XVG.write(("%d\t"*17) % tuple(self.frame_jtype))
        self.XVG.write("\n")
    def output_close(self):
        self.XVG.close()
        self.CXVG.close()
        if opts.verbose:
            sys.stderr.write("TOTAL:\t")
            sys.stderr.write(("%d\t"*4) % tuple(self.total_tally))
            if opts.qc:
                sys.stderr.write(("%d\t"*17) % tuple(self.total_jtype))
            sys.stderr.write("\n")


# Helper functions ############################################################
###############################################################################

def readndx(grps, ndx):
    readheaders = 0
    store_ats = False
    header = ''
    atoms = dict(zip(grps, [numpy.array([],int) for i in range(len(grps))]))
    ndxheader = re.compile('\s*\[\s*(\S(.*\S)?)\s*\]',re.IGNORECASE)
    for line in open(ndx).readlines():
        try:
            header = ndxheader.match(line).groups()[0].lower()
        except AttributeError:
            if store_ats:
                atoms[header] = numpy.append(atoms[header],map(int,line.split()))
            continue
        if header in atoms:
            store_ats = True
            readheaders += 1
            continue
        else:
            store_ats = False
            if readheaders == len(grps):
                break
    for k, ats in atoms.items():
        atoms[k] = atoms[k]-1
    return atoms

def sign(val):
    return int(val > 0) - int(val < 0)

#
##################################################################
##################################################################

if __name__ == "__main__":
    main(sys.argv)

