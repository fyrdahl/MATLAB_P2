
%number of dimensions
N = 2;
% Set directions
u = eye(N);

%Save starting position P0
P_0 = [ ];

%move to the minimum along direction
for i = 1:N
    P(i) = P_0 + u(i)
    
end

