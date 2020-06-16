function [sx sy] = defineStartingPoints
% define starting points for the streamlines

% cavity flow
crux = [0:0.1:1 , 0.5*ones(1,11); 0.5*ones(1,11), 0:0.1:1];
diagonals = [0:0.1:1,1:-0.1:0; 0:0.1:1, 0:0.1:1];
corner_dx = [ 0.9:0.03:1 ;...
    0.1:-0.03:0 ];
corner_sx = [0.1*ones(1,6) , 0:0.02:0.1, 0:0.02:0.1 ;...
    0:0.02:0.1 , 0:0.02:0.1 , 0.1*ones(1,6)];

corner_sx_up = [0.1*ones(1,6) , 0.9:0.02:1, 0:0.02:0.1 ;...
    0.9:0.02:1 , 0.9:0.02:1 , 0.9*ones(1,6)];

segm1 = [0:0.01:0.05; ...
    0.87*ones(1,6)];

segm2 = [0:0.005:0.05; ...
    0.80*ones(1,11)];

segm3 = [0:0.005:0.03; ...
    0.73*ones(1,7)];

s = [crux , diagonals , corner_dx , corner_sx, corner_sx_up, segm1, segm2, segm3];
sx = s(1,:);
sy = s(2,:);