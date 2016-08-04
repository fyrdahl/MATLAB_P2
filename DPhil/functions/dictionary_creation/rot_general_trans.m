function Rot_general = rot_general_trans(phi, theta)
% ROT_GENERAL_TRANS counter-clockwise rotation (degrees) about the y-axis
%
% Function to design a 3x3 rotation matrix, to rotate a given angle (degrees) about x, y and z.
%
% [ROT_GENERAL] = ROT_GENERAL_TRANS(PHI, THETA)
%
%   PHI is a scalar. The desired x,y rotation angle in degrees.
%   
%   THETA is a scalar. The desired z rotation angle in degrees.
%
%   ROT_GENERAL is a 3x3 matrix.
%
% Author: Jack Allen <jack.allen@jesus.ox.ac.uk>
% Copyright © 2016 University of Oxford

phi = deg2rad(phi);
theta = deg2rad(theta);

Rot_x(1,1:3) = ([1, 0 , 0]);
Rot_x(2,1:3) = ([0, cos(phi), -sin(phi)]);
Rot_x(3,1:3) = ([0, sin(phi) , cos(phi)]);

Rot_z1(1,1:3) = ([cos(theta), -sin(theta), 0 ]);
Rot_z1(2,1:3) = ([sin(theta) , cos(theta) , 0]);
Rot_z1(3,1:3) = ([0,0, 1]);

Rot_z2(1,1:3) = ([cos(-theta), -sin(-theta), 0 ]);
Rot_z2(2,1:3) = ([sin(-theta) , cos(-theta) , 0]);
Rot_z2(3,1:3) = ([0,0, 1]);

Rot_general = Rot_z1*Rot_x*Rot_z2;

end
