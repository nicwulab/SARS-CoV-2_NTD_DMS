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

hide all
show cartoon, spike
remove (hydro)
set cartoon_flat_sheets, 0
show sphere, spike and chain A and resi 232 and (not name c+o+n)
show sticks, spike and chain C and resi 355 and (not name c+o+n)
show sticks, spike and chain C and resi 466 and (not name c+o+n)
util.cnc all

set_view (\
     0.885324895,    0.424220085,   -0.190330520,\
     0.451583087,   -0.687021494,    0.569264531,\
     0.110734083,   -0.589938581,   -0.799813509,\
    -0.003420321,    0.001302402,  -39.754032135,\
   223.672531128,  243.811019897,  174.550552368,\
    26.438667297,   52.990764618,  -20.000000000 )
ray 1200, 1200
png ~\Desktop\G232.png
