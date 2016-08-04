function rot_z = rot_z(flipAngle)
% ROT_Z counter-clockwise rotation (degrees) about the z-axis
%
% Function to design a 3x3 rotation matrix, to rotate a given angle (degrees) about the z-axis.
%
% [ROT_Z] = ROT_Z(FA)
%
%   FA is a scalar. The desired rotation angle in degrees.
%
%   ROT_Z is a 3x3 matrix.
%
% Author: Jack Allen <jack.allen@jesus.ox.ac.uk>
% Copyright © 2016 University of Oxford
flipAngle = deg2rad(flipAngle);
rot_z = zeros(3,3);
rot_z(1,1:3) = ([cos(flipAngle), sin(flipAngle), 0 ]);
rot_z(2,1:3) = ([-sin(flipAngle) , cos(flipAngle) , 0]);
rot_z(3,1:3) = ([0,0, 1]);
end
