%  cap_pict.m
%    m file to generate coherency map pictures for 
%    CAP article

load c_map
m = [-20 30];  % orientation for 3D plot
mesh(20*c_map)
xlabel('Machine number')
ylabel('Machine number')
zlabel('Coherency measure')
view(m)
keyboard

load mode_map 
mesh(mode_map)
xlabel('x grid')
ylabel('y grid')
zlabel('Mode shape')
keyboard

