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
     0.885324895,    0.424220085,   -0.190330520,\
     0.451583087,   -0.687021494,    0.569264531,\
     0.110734083,   -0.589938581,   -0.799813509,\
    -0.003420321,    0.001302402,  -39.754032135,\
   223.672531128,  243.811019897,  174.550552368,\
    26.438667297,   52.990764618,  -20.000000000 )
    
ray 1200,1000
png ~/PDB/Images/spike_G232E.png, dpi=1200
