function Rot_general = Rot_general(phi, theta)
%clockwise rotation about the y-axis
phi = deg2rad(phi);
theta = deg2rad(theta);

Rot_x = zeros(3,3);
Rot_z = zeros(3,3);
Rot_general = zeros(3,3);

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
