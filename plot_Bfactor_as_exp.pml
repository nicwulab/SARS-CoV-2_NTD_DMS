set cartoon_transparency, 0
set ray_opaque_background, 1
set ray_shadows, 0
set valence, 0
bg_color white

load PDB/7b62_exp.pdb
create NTD, 7b62_exp and chain A and resi 14-301+1301
spectrum b, blue_white_red, minimum=0, maximum=100
color grey50, NTD and resi 14+15+22+30+31+36+39+45+46+47+54+55+61+62+63+70+71+76+77+78+85+86+93+94+101+102+103+104+105+109+110+117+125+133+134+141+154+155+173+189+205+214+260+290

load PDB/6zge.pdb
align 6zge and chain A, NTD and chain A
create spike, 6zge and (not (chain A and resi 14-301))
color grey50, spike

hide all
show surface, spike
color orange, NTD and resi 1301
util.cnc NTD and resi 1301
show cartoon, NTD

hide all
create RBD, 6zge and (resi 331-527)
create S2, 6zge and (not (resi 14-527))
create other_NTD, 6zge and chain A and chain B and resi 14-301
create other_NTD, 6zge and chain B+C and resi 14-301
color lightpink, other_NTD
color wheat, RBD
color palegreen, S2
show cartoon, RBD
show cartoon, NTD
show cartoon, other_NTD
show cartoon, S2

set_view (\
     0.796642780,    0.592633903,    0.118820176,\
     0.306436151,   -0.565446794,    0.765734971,\
     0.520990014,   -0.573608220,   -0.632077157,\
     0.000000000,    0.000000000, -218.880157471,\
    30.626884460,   -0.948181152,  -29.811660767,\
  -15204.945312500, 15642.705078125,  -20.000000000 )
ray 1000,1000
png PDB/Overlay_closeup_topdown.png