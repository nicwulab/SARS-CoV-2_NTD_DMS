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
show sticks, spike and chain C and resi 353 and (not name c+o+n)
show sticks, spike and chain C and resi 466 and (not name c+o+n)
util.cnc all

set_view (\
     0.721250474,    0.389912397,   -0.572497666,\
     0.650433660,   -0.097071558,    0.753327131,\
     0.238160446,   -0.915716708,   -0.323628455,\
    -0.004663851,    0.000188618,  -37.941143036,\
   225.039016724,  242.703659058,  172.866088867,\
    18.997615814,   56.717945099,  -20.000000000 )
ray 1200, 1200
png ~\Desktop\G232.png
