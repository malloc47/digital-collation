function final_matches = regressMatches(page1feat,page2feat, matches, thresh)
% REGRESSMATCHES     Score matches by what effect the match has on the
%                    resuting thin plate spline matrix transform, removing
%                    a fixed number of matches that score poorly.

%% Go through the usual steps to construct the thin plate spline matrix
% base_landmarks = page1feat;
% unregistered_landmarks = page2feat;

n=length(matches);
lambda = 1000;

% base_landmarks = fliplr(base_landmarks);
% unregistered_landmarks=fliplr(unregistered_landmarks);

vx = page1feat(:,1);
vy = page1feat(:,2);
ux = page2feat(:,1);
uy = page2feat(:,2);

L = computeWl(ux,uy);

L = pinv(L);
L = L(1:n,1:n);
Linv = pinv(lambda*L+eye(size(L,2),size(L,1)));

vhatx = Linv*vx;
vhaty = Linv*vy;

% diffx = vx - ux;
% diffy = vy - uy;

%% Compute bending energy and threshold

bending_energy = abs(abs(ux-vx) - abs(vx-vhatx)) + abs(abs(uy-vy) - abs(vy-vhaty));

bending_energy_sorted = sort((bending_energy),'ascend');

% index = floor(n*thresh); % floor(size(unregistered_landmarks,1)*tol);
index = floor(n-thresh); % Ignore thresh for right now.

% Select a score from the selection of scores to be the threshold
thresh_val = bending_energy_sorted(index);

final_matches = [];

% Only keep things below the selected threshold value (same as removing a
% fixed number of matches).
for i = 1:n
    if (bending_energy(i)) < thresh_val
        final_matches = [final_matches i];
    end
end

disp(['Regression removed ' num2str(size(matches,2)-size(final_matches,2)) ]);

end

% % Unused
% function normed = normvect(vect1,vect2)
% 
% normed = (vect1(:,1).^2 + vect2(:,1).^2).^0.5 + (vect1(:,2).^2 + vect2(:,2).^2).^0.5;
% 
% end


function [wL]=computeWl(xp, yp)

np = length(xp);

rXp = repmat(xp(:),1,np); % 1xNp to NpxNp
rYp = repmat(yp(:),1,np); % 1xNp to NpxNp

wR = sqrt((rXp-rXp').^2 + (rYp-rYp').^2); % compute r(i,j)

wK = radialBasis(wR); % compute [K] with elements U(r)=r^2 * log (r^2)
wP = [ones(np,1) xp(:) yp(:)]; % [P] = [1 xp' yp'] where (xp',yp') are n landmark points (nx2)
wL = [wK wP;wP' zeros(3,3)]; % [L] = [[K P];[P' 0]]

end

function [ko]=radialBasis(ri)

r1i = ri;
r1i(find(ri==0))=realmin; % Avoid log(0)=inf
ko = 2*(ri.^2).*log(r1i);

end
