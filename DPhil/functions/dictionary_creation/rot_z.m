function Rot_z = Rot_z(flipAngle)
%produces a clockwise rotation matrix, about the z  
flipAngle = deg2rad(flipAngle);
Rot_z = zeros(3,3);
Rot_z(1,1:3) = ([cos(flipAngle), -sin(flipAngle), 0 ]);
Rot_z(2,1:3) = ([sin(flipAngle) , cos(flipAngle) , 0]);
Rot_z(3,1:3) = ([0,0, 1]);
end
