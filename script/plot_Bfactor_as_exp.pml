set cartoon_transparency, 0
set ray_opaque_background, 1
set ray_shadows, 0
set valence, 0
bg_color white

load PDB/7b62_exp.pdb
create NTD, 7b62_exp and chain A and resi 14-301+1301
spectrum b, blue_white_red, minimum=0, maximum=100
color grey50, NTD and resi 14+15+22+23+30+31+39+45+46+47+54+55+61+62+63+70+71+76+77+78+85+93+94+96+101+102+104+109+117+125+133+134+141+154+155+173+187+205+214+260+290

load PDB/6zge.pdb
align 6zge and chain A, NTD and chain A
create spike, 6zge and (not (chain A and resi 14-301))
color grey50, spike

hide all
show surface, spike
color orange, NTD and resi 1301
util.cnc NTD and resi 1301
show cartoon, NTD
