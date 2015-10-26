function Rot_y = Rot_y(flipAngle)
%clockwise rotation about the y-axis
flipAngle = deg2rad(flipAngle);
Rot_y = zeros(3,3);
Rot_y(1,1:3) = ([cos(flipAngle), 0 , sin(flipAngle)]);
Rot_y(2,1:3) = ([0 , 1 , 0]);
Rot_y(3,1:3) = ([-sin(flipAngle), 0 , cos(flipAngle)]);

end
