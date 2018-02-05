# fluxer
A script for counting crossing events across a finite width xy-plane membrane.

The system is assumed periodic in z but the pbc in the trajectory must have been
removed (trjconv option -pbc nojump).
Two index groups define the top and bottom of the membrane. A crossing is only
counted when an atom exits on the opposite side of the membrane to the one it
came in.
When analyzing a single channel, the option -mult can be passed to restrict
counting to the cylinder of radius mult\*RoG, where RoG is the radius of gyration
(not mass weighted) of the delimiting groups. Care must be taken to ensure these
groups remain whole for the whole trajectory.
fluxer.py assumes an xy flat membrane slab. Other topologies or heavy
fluctuations will produce erroneous results.

The schematic of the counting is as follows:
```
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
```
An event is recorded if a water molecule arrives at the water slab of a box,
having last been at the water slab of an adjacent box. In addition, the crossing
must be in line with the channel, meaning the molecule must have been in the
channel ('eligible') in the frame preceding the one when it reached the
water slab.

Columns in the output .xvg and the analysis total are as follows:
```
1: Time
2: +Flux
3: -Flux
4: Crossings that were bigjumps
5: Total bigjumps
6-22: Jumptype discrimination.
```
Bigjump signals crossings between non-consecutive regions (namely, directly
jumping from one water slab to the next in one step). Having a large number of
these may indicate too large a frame step for the analysis. 300 ps per frame
is a typically safe value (for water).

The jumptypes code for different events that are recorded. These are numbered
0-16 and printed out in order in the .xvg columns 6-22.
Each frame the slab and eligibility of each molecule is compared to the previous
frame's. The jumptypes codes are as per the following matrix (the 'e E n N'
labels are as per the legend above):
```
         To(new):   
From:     e   E   n   N    +: stands for a clean crossing;                
(old) e *16 +15 *14 +13    *: is a crossing that can only happen
      E ?12 *11 ?10  *9       together with a bigjump;
      n  x8  x7  x6  x5    x: will not be counted as a crossing.
      N  ?4  *3  x2  x1    ?: is a crossing that should never be observed
                              (out of internal consistency - it requires
                              a jump so big that it is prevented by the
                              pbc treatment of nojump).
```

In addition to these, jumptype = 0 means no slab difference, or a membrane-water
transition back to the last water slab the molecule visited.
Eligibility is defined as being within the z-axis-oriented cylinder defined by
the ref groups' CoG and RoG (scaled by the -mult option). 
The assignment of `*` and `x` is somewhat arbitrary regarding some exotic jumps
that will only be of concern when the frame step is so high that the results
will be rubbish anyway. The bottom line is that one should aim for mostly
`+` events; too many bigjumps indicate that the step is inappropriate.
