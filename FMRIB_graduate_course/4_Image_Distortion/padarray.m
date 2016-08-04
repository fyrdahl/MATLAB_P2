function obj_padded = padarray( object, pad )

%Replace padarray function
%Oct 2013
%Mark Chiew

obj_padded = zeros(size(object)+2*pad);
obj_padded(pad(1)+1:pad(1)+size(object,1),pad(2)+1:pad(2)+size(object,2)) = object;

end

