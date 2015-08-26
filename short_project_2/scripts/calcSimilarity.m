%% simularity measure
for n = 1:2
    % Cosine similarity (does not depend on magnitude of each vector)
    similarity(n,T1,T2) = dot(ySim,y(n,:))/(norm(ySim)*norm(y(n,:)));
end
%     end
% end