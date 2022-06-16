set ray_shadows, 0
bg_color white
set specular_intensity, 0.2
set valence, 0
set seq_view, 1

load PDB/S50Q_fastrelax.pdb, spike

color gray90, spike and chain A
color lightblue, spike and chain B
color pink, spike and chain C

hide all
show cartoon, spike
remove (hydro)
set cartoon_flat_sheets, 0

show sticks, spike and chain A and resi 37 and (not name c+o+n)
show sticks, spike and chain A and resi 256 and (not name c+o+n)
show sticks, spike and chain B and resi 1811 and (not name c+o+n)
util.cnc all

distance dist1, /spike//A/GLN`37/NE2, /spike//A/THR`256/OG1
distance dist1, /spike//A/GLN`37/OE1, /spike//B/SER`1811/OG
hide label, dist1
color black, dist1

set_view (\
    -0.057936318,   -0.984386742,   -0.166202441,\
    -0.996987939,    0.065635577,   -0.041215103,\
     0.051479261,    0.163313538,   -0.985224605,\
     0.001402628,   -0.000129504,  -56.949554443,\
   225.841873169,  236.471176147,  210.875091553,\
    35.989475250,   78.681236267,  -20.000000000 )

ray 1200,1000
png ~/PDB/Images/spike_S50Q.png, dpi=1200