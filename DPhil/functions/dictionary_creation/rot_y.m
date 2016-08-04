<<<<<<< HEAD
function rot_y = rot_y(FA)
% ROT_Y counter-clockwise rotation (degrees) about the y-axis
%
% Function to design a 3x3 rotation matrix, to rotate a given angle (degrees) about the y-axis.
%
% [ROT_Y] = ROT_Y(FA)
%
%   FA is a scalar. The desired rotation angle in degrees.
%
%   ROT_Y is a 3x3 matrix.
%
% Author: Jack Allen <jack.allen@jesus.ox.ac.uk>
% Copyright © 2016 University of Oxford

rot_y = zeros(3,3);
rot_y(1,1:3) = ([cosd(FA), 0 , -sind(FA)]);
rot_y(2,1:3) = ([0 , 1 , 0]);
rot_y(3,1:3) = ([sind(FA), 0 , cosd(FA)]);
=======
function Rot_y = Rot_y(flipAngle)
%clockwise rotation about the y-axis
flipAngle = deg2rad(flipAngle);
Rot_y = zeros(3,3);
Rot_y(1,1:3) = ([cos(flipAngle), 0 , sin(flipAngle)]);
Rot_y(2,1:3) = ([0 , 1 , 0]);
Rot_y(3,1:3) = ([-sin(flipAngle), 0 , cos(flipAngle)]);
>>>>>>> 47303176e9617d0556d5fc1882e1cda83864dfe1

end
