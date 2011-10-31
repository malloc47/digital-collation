function [tform,K]=findTransform(landmarks1, landmarks2)
% FINDTRANSFORM     Returns the thin plate transform derived from the
%                   provided landmarks.

[m,n] = size(landmarks1);
[m1,n1] = size(landmarks2);

% Error cases
if(m ~= m1 || n ~= n1)
    disp('Warning: Mismatch in number of landmarks');
    return;
end
if(n ~= 2)
    disp('Warning: Landmarks are not 2D');
    return;
end

% Make tps matrix from landmarks1
L = zeros(m,m);
for i=1:m
    for j=1:m
        d = landmarks1(i,:) - landmarks1(j,:);
        r = d * d';
        r = sqrt(r);
        if(r > 0.0)
            L(i,j) = r*r*log(r);
        end
        
    end
end

L = L + L';
L = L ./ 2;

K = [L ones(m,1) landmarks1; ([ones(m,1) landmarks1])' zeros(n+1,n+1) ];

% Find transform factors using landmarks2
tempL2 = [landmarks2; zeros(n+1,n)];
tform = pinv(K)*tempL2;


disp(['Translation: ' num2str(tform(m+1,:))]);

disp(['Affine: ' num2str(tform(m+2,:)) ' ; ' num2str(tform(m+3,:))]);

end