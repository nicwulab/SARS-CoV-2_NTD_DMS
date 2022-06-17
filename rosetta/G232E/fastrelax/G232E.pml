set ray_shadows, 0
bg_color white
set specular_intensity, 0.2
set valence, 0
set seq_view, 1

load PDB/G232E_fastrelax.pdb, spike

color gray90, spike and chain A
color lightblue, spike and chain B
color pink, spike and chain C

hide all
show cartoon, spike
remove (hydro)
set cartoon_flat_sheets, 0

show sticks, spike and chain A and resi 214 and (not name c+o+n)
show sticks, spike and chain C and resi 2531 and (not name c+o+n)
show sticks, spike and chain C and resi 2644 and (not name c+o+n)
util.cnc all

distance dist1, /spike//A/GLU`214/OE1, /spike//C/ARG`2644/NE
distance dist1, /spike//A/GLU`214/OE1, /spike//C/ARG`2644/NH1
distance dist1, /spike//A/GLU`214/OE2, /spike//C/TRP`2531/NE1
hide label, dist1
color black, dist1

set_view (\
     0.847988725,    0.231995121,   -0.476535380,\
     0.493765056,   -0.019048700,    0.869381487,\
     0.192616314,   -0.972525597,   -0.130706027,\
     0.000000000,    0.000000000,  -46.132404327,\
   225.516006470,  244.415618896,  171.548156738,\
    16.630271912,   75.634521484,  -20.000000000 )
    
ray 1200,1000
png ~/PDB/Images/spike_G232E.png, dpi=1200
