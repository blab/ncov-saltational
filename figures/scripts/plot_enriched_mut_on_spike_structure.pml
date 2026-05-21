set assembly, 1
fetch 7KRQ
hide everything, 7KRQ
show cartoon, 7KRQ
color 0xd3d3d3
set seq_view, 1
bg_color white
select enriched_pos_RBD, 7KRQ and resi 332+339+340+346+356+368+373+375+376+385+403+417+440+444+445+446+449+450+452+460+477+478+484+486+490+493+494+496+498+500+501+505 and name CA
show spheres, enriched_pos_RBD
color 0x00cdcd, enriched_pos_RBD
select enriched_pos_NTD, 7KRQ and resi 19+64+83+95+136+152+157+164+210+213+214+252+257 and name CA
show spheres, enriched_pos_NTD
color 0xcd1076, enriched_pos_NTD
select enriched_pos_other_S1, 7KRQ and resi 568+655+681 and name CA
show spheres, enriched_pos_other_S1
color 0x483d8b, enriched_pos_other_S1
select enriched_pos_S2, 7KRQ and resi 764+936+954+981+1002+1118+1190 and name CA
show spheres, enriched_pos_S2,
color 0xee9a00, enriched_pos_S2
select enriched_pos_other, 7KRQ and resi 9 and name CA
show spheres, enriched_pos_other
color 0x000000, enriched_pos_other
set sphere_scale, 1.5

orient
rotate z, -90
rotate y, 30
ray 2400,2400