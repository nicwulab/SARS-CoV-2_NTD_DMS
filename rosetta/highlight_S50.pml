fetch 6zge
set_name 6zge,spike
set ray_shadows, 0
bg_color white
set specular_intensity, 0.2
set valence, 0
set seq_view, 1

color gray90, spike and chain A
color lightblue, spike and chain B
color pink, spike and chain C
cartoon loop, resi 40-60
extract S2, spike and chain B+C and resi 532-1000

hide all
show cartoon, spike
remove (hydro)
set cartoon_flat_sheets, 0
show surface, S2
show sticks, chain A and resi 50+304 and (not name c+n+o)
util.cnc all

color lightblue, S2
distance dist1, /spike//A/LYS`304/NZ, /spike//A/SER`50/OG
color black, dist1
hide label, dist1

set_view (\
     0.010553611,   -0.546536744,   -0.837353528,\
    -0.902329862,   -0.366049409,    0.227549508,\
    -0.430876851,    0.753180325,   -0.497027338,\
     0.011521675,    0.001096303,  -31.734323502,\
   223.531661987,  240.716476440,  211.528900146,\
    19.745407104,   49.507717133,  -20.000000000 )

ray 1200, 1200
png ~\Desktop\S50.png

