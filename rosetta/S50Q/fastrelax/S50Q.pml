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
-0.094882026,   -0.979946375,   -0.175215155,\
-0.930556595,    0.149831653,   -0.334074795,\
0.353626907,    0.131349936,   -0.926112652,\
0.000000000,    0.000000000,  -44.250446320,\
223.780212402,  238.443832397,  206.480911255,\
27.898347855,   60.602527618,  -20.000000000 )

ray 1200,1000
png ~/PDB/Images/spike_S50Q.png, dpi=1200
