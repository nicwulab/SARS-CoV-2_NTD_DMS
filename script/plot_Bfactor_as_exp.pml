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

set ray_trace_mode,  0
set_view (\
     0.619438767,   -0.098452665,    0.778835654,\
    -0.594106674,    0.589699030,    0.547061801,\
    -0.513138950,   -0.801594436,    0.306791097,\
    -0.001331814,    0.000099331, -555.848693848,\
    36.505077362,  -48.433963776,  -21.978576660,\
  -13867.971679688, 14979.678710938,  -20.000000000 )
ray 1200,1200
png PDB/Overlay_front.png

set cartoon_transparency, 0.5, RBD
set cartoon_transparency, 0.5, other_NTD
set cartoon_transparency, 0.5, S2
set ray_trace_mode,  1

set_view (\
     0.501908720,    0.052816618,    0.863294184,\
    -0.639454186,    0.694736421,    0.329270452,\
    -0.582373142,   -0.717311084,    0.382471681,\
    -0.000407148,   -0.000369415, -244.470108032,\
    13.404919624,    7.596892357,  -22.284166336,\
  -14179.359375000, 14668.291015625,  -20.000000000 )
ray 1000,1000
png PDB/Overlay_closeup_0.png


set_view (\
     0.658307672,    0.273608267,   -0.701245010,\
     0.441339970,    0.614371240,    0.654023051,\
     0.609776318,   -0.740047932,    0.283686966,\
    -0.000407148,   -0.000369415, -218.851196289,\
    13.404919624,    7.596892357,  -22.284166336,\
  -14204.978515625, 14642.671875000,  -20.000000000 )
ray 1000,1000
png PDB/Overlay_closeup_180.png

set_view (\
     0.796642780,    0.592633903,    0.118820176,\
     0.306436151,   -0.565446794,    0.765734971,\
     0.520990014,   -0.573608220,   -0.632077157,\
     0.000000000,    0.000000000, -218.880157471,\
    30.626884460,   -0.948181152,  -29.811660767,\
  -15204.945312500, 15642.705078125,  -20.000000000 )
ray 1000,1000
png PDB/Overlay_closeup_topdown.png
