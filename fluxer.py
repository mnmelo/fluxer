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
import numpy as np
import mdreader

# Static descriptions ##################################################
########################################################################

desc = '''
Counts crossing events across a finite width xy-plane membrane.
The system is assumed periodic in z but the pbc in the trajectory must have
been removed (trjconv option -pbc nojump).
Two index groups define the top and bottom of the membrane. A crossing is only
counted when an atom exits on the opposite side of the membrane to the one it
came in.
When analyzing a single channel, the option -mult can be passed to restrict
counting to the cylinder of radius mult*RoG, where RoG is the radius of
gyration (not mass weighted) of the delimiting groups. Care must be taken to
ensure these groups remain whole for the whole trajectory.
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
having last been at the water slab of an adjacent box. In addition, the
crossing must be in line with the channel, meaning the molecule must have been
in the channel ('eligible') in the frame preceding the one when it reached the
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
Each frame the slab and eligibility of each molecule is compared to the
previous frame's. The jumptypes codes are as per the following matrix (the
'e E n N' labels are as per the legend above):

         To(new):   
From:     e   E   n   N    +: stands for a clean crossing;                
(old) e *16 +15 *14 +13    *: is a crossing that can only happen
      E ?12 *11 ?10  *9       together with a bigjump;
      n  x8  x7  x6  x5    x: will not be counted as a crossing.
      N  ?4  *3  x2  x1    ?: is a crossing that should never be observed
                              (out of internal consistency - it requires
                              a jump so big that it is prevented by the
                              pbc treatment of nojump).

In addition to these, jumptype = 0 means no slab difference, or a
membrane-water transition back to the last water slab the molecule visited.
Eligibility is defined as being within the z-axis-oriented cylinder defined by
the ref groups' CoG and RoG (scaled by the -mult option). 
The assignment of '*' and 'x' is somewhat arbitrary regarding some exotic jumps
that will only be of concern when the frame step is so high that the results
will be rubbish anyway. The bottom line is that one should aim for mostly
'+' events; too many bigjumps indicate that the step is inappropriate.
   
'''

progname = "fluxer.py"
vernum = "0.2"

ver = """%s v%s
by Manuel Melo (m.n.melo@itqb.unl.pt)""" % (progname, vernum)

header = '''#   Produced by fluxer.py %s (written by Manuel Melo; m.n.melo@itqb.unl.pt)
#   %s
#\n''' #to be formatted with % (vernum, " ".join(sys.argv))

fluxheader   = "#\n#   Flux vs time.\n"
fluxprologue = "# Columns: 1:Time, 2:Upward flux, 3:Downward flux, 4:Bigjumps counted as flux, 5:Total bigjumps, 6-22:Jumptype discrimination.\n"

cumfluxheader   = "#\n#   Cumulative flux vs time.\n"
cumfluxprologue = "# Columns: 1:Time, 2:Total cum. flux, 3:Upward cum. flux, 4:Downward cum. flux.\n"

def main():
    topring, botring = syst.ndxgs
    if not len(topring):
        ValueError("Empty or not-found index top group \'%s\'."
                   % syst.opts.topgrp)
    if not (len(botring) or opts.toroidal):
        ValueError("Empty or not-found index bottom group \'%s\'."
                   % syst.opts.botgrp)
    bothrings = topring + botring
    ringatms_inv = 1./len(bothrings)
    waters = None
    if syst.opts.toroidal and syst.opts.tor_scale:
        rel_membrane_top = syst.opts.toroidal/(2*syst.dimensions[2])
        rel_membrane_bot = -rel_membrane_top 

    waters = WaterCollection(syst)

    for frame in syst.iterate(p=1): #must be serial

        x_dim, y_dim, z_dim = frame.dimensions[:3]
        if syst.opts.toroidal:
            channel_cog = topring.center_of_geometry()
            if not syst.opts.tor_scale:
                rel_membrane_top = opts.toroidal/(2*self.z_dim)
                rel_membrane_bot = -rel_membrane_top 
        else:
            top_cog = topring.center_of_geometry()
            bot_cog = botring.center_of_geometry()
            channel_cog = (top_cog + bot_cog)/2
            # xy-plane radius of gyration; not mass weighted
            cutoffsq = syst.opts.radius_multsq * ringatms_inv * np.sum(
                (bothrings.positions[:,:2] - channel_cog[:2])**2)
            rel_membrane_top = (top_cog[2] - channel_cog[2])/z_dim
            rel_membrane_bot = -rel_membrane_top 

        if not syst.iterframe:
            # The first frame is an initializer
            waters.update_flux(frame, channel_cog, rel_membrane_top,
                               cutoffsq, init=True)
            continue
        else:
            waters.update_flux(frame, channel_cog, rel_membrane_top, cutoffsq)
        waters.output_frame(frame)

    waters.output_close()

# Classes ##############################################################
########################################################################

class WaterCollection:
#    #                                     Counted  Total    
#    #                         +Flux -Flux Bigjumps Bigjumps  
#    frame_tally = numpy.array([  0,    0,     0,      0   ])

    def __init__(self, syst):
        self.syst = syst
        self.water_atms = syst.select_atoms("name %s"% syst.opts.watername)
        if not self.water_atms:
            ValueError("No atoms found with name '%s'." % syst.opts.watername)
        nwaters = len(self.water_atms)

        self.XVG = open(syst.opts.outfile, 'w')
        self.CXVG = open(syst.opts.cumfile, 'w')
        if syst.opts.xvgheaders:
            self.XVG.write(fluxheader)
            self.XVG.write(header % (vernum, " ".join(syst.arguments)))
            self.XVG.write(fluxprologue)
            self.CXVG.write(cumfluxheader)
            self.CXVG.write(header % (vernum, " ".join(syst.arguments)))
            self.CXVG.write(cumfluxprologue)

        self.total_tally = np.zeros(4, dtype=np.int) 
        self.total_jtype = np.zeros(17, dtype=np.int)
        self.jumptype_allow = np.array(
            [0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1])
        self.inmembrane = mdreader.ThenNow()
        self.slab = mdreader.ThenNow()
        self.waterslab = mdreader.ThenNow()
        self.waterslab.fill(np.zeros(nwaters, dtype=np.int))
        self.eligible = mdreader.ThenNow(True, True)
        self.inmembrane = mdreader.ThenNow()

    def update_flux(self, frame, channel_cog, rel_top, cutoffsq, init=False):
        # bring in the centered coordinates
        self.coords = self.water_atms.positions - channel_cog

        # set the slab and waterslab
        z_dim = frame.dimensions[2]
        rx = self.coords[:,2]/z_dim
        mdf = np.modf(rx)
        self.inmembrane.update(np.abs(mdf[0]-np.round(mdf[0])) < rel_top, init)
        self.slab.update((np.round(rx) * 2*self.inmembrane.new +
                          (mdf[1]*2+np.sign(rx)) * (~self.inmembrane.new))
                         .astype(np.int), init)
        slabdiff = self.slab.new - self.slab.old
        bigjump = np.abs(slabdiff) > 1
        wslab_updt0 = ~bigjump * self.inmembrane.new
        wslab_updt1 = bigjump * self.inmembrane.new
        self.waterslab.update(self.waterslab.new * wslab_updt0 +
                              self.slab.new * (~self.inmembrane.new) +
                              (self.slab.new - np.sign(slabdiff)) *
                              wslab_updt1,
                              init)
        wslabdiff = self.waterslab.new - self.waterslab.old

        # set eligibility
        if self.syst.opts.radius_multsq:
            xy = self.coords[:,:2] - np.round(
                self.coords[:,:2]/frame.dimensions[:2]) * frame.dimensions[:2]
            self.eligible.update((xy**2).sum(axis=1) < cutoffsq, init)

        # jumptype and flux
        jumptype =  ((self.waterslab.old != 0) * (wslabdiff != 0) *
                     (np.left_shift(self.eligible.old, 3) +
                      np.left_shift(self.inmembrane.old, 2) +
                      np.left_shift(self.eligible.new, 1) +
                      self.inmembrane.new + 1))
        flux = np.sign(wslabdiff) * self.jumptype_allow[jumptype]

        # tally up
        self.frame_jtype = np.bincount(jumptype, minlength = 17)
        #values:
        # ((+Flux, -Flux, Bigjumps counted as jumps, Total bigjumps), jumptype)
        #  jumptype reports jump type, as per the from:to matrix above.
        #  jumptype = 0 means no slab difference (not the same as flux = 0 as
        #  there may be a slab difference but no jump, as in a
        #  notelig+even -> notelig+noteven transition).
        self.frame_tally = np.array([np.sum(flux>0),
                                     np.sum(flux<0),
                                     np.sum(bigjump),
                                     np.sum(np.multiply(bigjump,flux!=0))])
        self.total_tally += self.frame_tally
        self.CXVG.write("%f\t%d\t%d\t%d\n" % (frame.time,
                                              self.total_tally[0] +
                                                self.total_tally[1],
                                              self.total_tally[0], 
                                              self.total_tally[1]))
        self.total_jtype += self.frame_jtype

        if self.syst.opts.verbose > 1:
            jump         = (self.waterslab.old!=0)*(wslabdiff!=0)
            # Add here water indices to be debugged every frame.
            for mol in []:
                jump[mol]    = True
            if np.sum(jump):
                sys.stderr.write("Frame: %d   %.1f ps  rel.top: %f  zdim: %f\n"
                                 % (frame.frame, frame.time, rel_top, z_dim))
            dbgindices   = np.where(jump)
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
   bigjump: %r  jumptype: %d\n"""
                    % (dbgindices[0][tgt_mol],
                       dbgcoords[tgt_mol],
                       dbginmemb_o[tgt_mol],
                       dbginmemb_n[tgt_mol],
                       dbgslab_o[tgt_mol],
                       dbgslab_n[tgt_mol],
                       dbgslab_d[tgt_mol],
                       dbgwslab_o[tgt_mol],
                       dbgwslab_n[tgt_mol],
                       dbgwslab_d[tgt_mol],
                       dbgbjump[tgt_mol],
                       dbgjtype[tgt_mol]))

    def output_frame(self, frame):
        self.XVG.write("%f\t" % (frame.time))
        self.XVG.write(("%d\t"*4) % tuple(self.frame_tally))
        if self.syst.opts.qc:
            self.XVG.write(("%d\t"*17) % tuple(self.frame_jtype))
        self.XVG.write("\n")

    def output_close(self):
        self.XVG.close()
        self.CXVG.close()
        if self.syst.opts.verbose:
            sys.stderr.write("TOTAL:\t")
            sys.stderr.write(("%d\t"*4) % tuple(self.total_tally))
            if self.syst.opts.qc:
                sys.stderr.write(("%d\t"*17) % tuple(self.total_jtype))
            sys.stderr.write("\n")


#
##################################################################
##################################################################

if __name__ == "__main__":
    syst = mdreader.MDreader(description=ver)
    ## Options
    #syst.add_argument('-o', metavar='OUT', dest='outfile', default='flux.xvg',
    #    help = 'file\tThe flux vs time file. Read the --help for the meaning '
    #                 'of each column in the output.')
    syst.add_argument('-oc', metavar='OUT', dest='cumfile',
        default='flux-cumul.xvg',
        help = 'file\tThe cumulative flux vs time file. Columns are: time, '
                     'total flux, upward flux, downward flux.')
    syst.add_argument('-atname', metavar='NAME', dest='watername', default='W',
        help = 'str\tThe atom name counted for the crossings.')
    syst.add_argument('-toroidal', metavar='THICKNESS', type=float,
        dest='toroidal', default = 0.,
        help = 'real\tAssumes a toroidal pore and sets its z-thickness, '
                    'ignoring the other top/bottom limits. Only the "top" '
                    'index group will be used and its CoM should then match '
                    'the membrane\'s (use all the lipids, for example). Such '
                    'channel-independent pores should always be analyzed with '
                    '-mult set to 0, as no pore tracking is implemented '
                    '(yet); no distinction will be made between '
                    'through-the-membrane and through-the-pore crossings.')
    syst.add_argument('-noscale', action='store_false', dest='tor_scale',
        help = 'bool\tBy default the thickness defined with -toroidal is '
                     'scaled with the box\'s Z dimension. Set to keep '
                     'constant.')
    syst.add_argument('-mult', metavar='MULT', type=lambda x:float(x)**2,
        dest='radius_multsq', default=0,
        help = 'real\tMultiplier of the xy radius of gyration of the top and '
                     'bottom groups to consider as part of the channel. Set '
                     'to 0 (default) to count all membrane crossings.')
    syst.add_argument('-noqc', action='store_false', dest='qc',
        help = 'bool\tBy default quality control is done by outputting '
                     'individual counts of crossing types. Set to disable.')
    syst.add_argument('-nohdr', action='store_false', dest='xvgheaders',
        help = 'bool\tVersion info, used command line parameters, and column '
                     'identification are added (#commented) to the header of '
                     'the output .xvg files. Set to disable header inclusion.')
    syst.add_ndx(ndxparms=['Top group (or whole membrane if using -toroidal)',
                           'Bottom group (ignored if using -toroidal)'])
    syst.setargs(o='flux.xvg')

    #mdreader.check_outfile(syst.opts.outfile) # already done by default
    mdreader.check_outfile(syst.opts.cumfile)
    mdreader.check_positive(syst.opts.toroidal)
    if len(syst.opts.watername) > 5 or ' ' in syst.opts.watername:
        raise ValueError('Invalid atom name %s' % opts.watername)
    if syst.opts.toroidal and syst.opts.radius_multsq:
        raise NotImplementedError("You are analyzing a toroidal pore and have "
                                  "set the -mult option. This is currently "
                                  "not implemented because the "
                                  "pore-delimiting atoms may change along the "
                                  "trajectory. Use -mult 0 (the default) when "
                                  "using -toroidal. The downside is that all "
                                  "crossings will be counted, even those "
                                  "outside the pore, which will hopefully be "
                                  "negligible.")

    # some more checking is done along the way.
    
    ## Let rip
    main()

