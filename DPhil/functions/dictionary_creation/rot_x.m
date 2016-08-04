<<<<<<< HEAD
function rot_x = rot_x(flipAngle)
% ROT_X counter-clockwise rotation (degrees) about the x-axis
%
% Function to design a 3x3 rotation matrix, to rotate a given angle (degrees) about the x-axis.
%
% [ROT_X] = ROT_X(FA)
%
%   FA is a scalar. The desired rotation angle in degrees.
%
%   ROT_X is a 3x3 matrix.
%
% Author: Jack Allen <jack.allen@jesus.ox.ac.uk>
% Copyright © 2016 University of Oxford
rot_x = zeros(3,3);
rot_x(1,1:3) = ([1, 0 , 0]);
rot_x(2,1:3) = ([0, cosd(flipAngle), sind(flipAngle)]);
rot_x(3,1:3) = ([0, -sind(flipAngle) , cosd(flipAngle)]);
=======
function Rot_x = Rot_x(flipAngle)
%clockwise rotation about the x-axis
flipAngle = deg2rad(flipAngle);
Rot_x = zeros(3,3);
Rot_x(1,1:3) = ([1, 0 , 0]);
Rot_x(2,1:3) = ([0, cos(flipAngle), -sin(flipAngle)]);
Rot_x(3,1:3) = ([0, sin(flipAngle) , cos(flipAngle)]);
>>>>>>> 47303176e9617d0556d5fc1882e1cda83864dfe1

end
