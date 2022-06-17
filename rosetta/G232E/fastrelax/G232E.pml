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
     0.766683638,    0.341676295,   -0.543545187,\
     0.600902498,   -0.083806351,    0.794909656,\
     0.226051256,   -0.936069310,   -0.269570023,\
    -0.004621752,    0.000350639,  -37.606369019,\
   225.725814819,  242.874969482,  172.056015015,\
    22.679824829,   52.480175018,  -20.000000000 )
    
ray 1200,1000
png ~/PDB/Images/spike_G232E.png, dpi=1200
