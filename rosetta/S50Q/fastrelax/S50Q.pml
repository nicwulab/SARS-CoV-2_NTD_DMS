set ray_shadows, 0
bg_color white
set specular_intensity, 0.2
set valence, 0
set seq_view, 1

load PDB/S50Q_fastrelax.pdb, spike

color gray90, spike and chain A
color lightblue, spike and chain B
color pink, spike and chain C
cartoon loop, resi 26-46
extract S2, spike and chain B+C and resi 532-1000

hide all
show cartoon, spike
remove (hydro)
set cartoon_flat_sheets, 0
show surface, S2
show sticks, chain A and resi 37+286 and (not name c+n+o)
util.cnc all

color lightblue, S2
distance dist1, /spike//A/LYS`286/NZ, /spike//A/GLN`37/OE1
color black, dist1
hide label, dist1


set_view (\
     0.063089155,   -0.537977219,   -0.840579927,\
    -0.915875852,   -0.365781099,    0.165365636,\
    -0.396430761,    0.759446681,   -0.515805066,\
     0.011393048,    0.001124245,  -30.851274490,\
   223.444351196,  237.225250244,  208.300979614,\
    23.945901871,   45.307220459,  -20.000000000 )

ray 1200,1000
png ~/PDB/Images/spike_S50Q.png, dpi=1200
