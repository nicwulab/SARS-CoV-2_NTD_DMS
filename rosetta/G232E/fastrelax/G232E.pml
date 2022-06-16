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
     0.786634743,    0.186685786,   -0.588513851,\
     0.600898981,   -0.012520080,    0.799221694,\
     0.141836345,   -0.982335687,   -0.122029826,\
    -0.002964417,    0.000841832,  -68.545455933,\
   241.684661865,  223.567001343,  175.212173462,\
    26.075973511,  110.773468018,  -20.000000000 )
    
ray 1200,1000
png ~/PDB/Images/spike_G232E.png, dpi=1200