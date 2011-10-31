function final_overlay = constructOverlay(registered,base,alpha)
% CONSTRUCTOVERLAY     Creates a composite image from the two input images
%                      along with an optional alpha value to change the
%                      transparency.

if nargin < 3
    alpha = 0.5;
end

registered = imresize(registered,[size(base,1) size(base,2)]);

final_overlay = base;
final_overlay(:,:,1) = alpha.*base(:,:,1) + (1-alpha).*registered(:,:,1);
final_overlay(:,:,2) = alpha.*base(:,:,2) + (1-alpha).*registered(:,:,2);
final_overlay(:,:,3) = alpha.*base(:,:,3) + (1-alpha).*registered(:,:,3);

end