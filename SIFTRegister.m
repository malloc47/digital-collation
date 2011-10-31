function registered = SIFTRegister(page1file, page2file,param)
%SIFTRegister   Register template (page2file) to target (page1file)
%   registered = SIFTRegister(page1file,page2file,param) returns an image 
%                of the registered template.
%   page1file -> The file name (absolute or relative) to the image to be registered
%   page2file -> The file name (absolute or relative) to the image that 
%                page1file will be registered to
%   param is a cell array of the following form:
%   param{1} = SIFT magnification (default: 8)
%   param{2} = Start threshold for matching (default: 10)
%   param{3} = Target number of total matches (default: 500)
%   param{4} = Number of matches that regression will remove (default: 3)
%   param{5} = [x,y] where x,y are the number of slices the page will be
%   divided into for doing the neighborhood consistency check. (default: [3,6])
%   param{6} = % of low-scoring matches to remove from each cell (default: 0.65)
%   param{7} = Target number of matches to randomly sample down to, if the
%   total number of matches is greater than this target (default: 0, disabled)
%   param{8} = Boolean to toggle display of figures (default: false)
%   param{9} = Boolean to toggle interface to adjust matches (may require a
%              long time to initialize interface if given many matches)(default: false)
%   
%   This process can be broken into several discrete phases:
%   1. Initialization and extraction of SIFT points
%   2. Iterative matching of sift points until a target number of matches is achieved
%   3. Removal of poor matches using multiple techniques
%   4. Registration of pages using thin plate spline warping

%   Copyright � 2009-2010 Center for Digital Humanities and the Department of 
%   Computer Science at the University of South Carolina

% Default parameters if none are provided
if nargin < 3
    %        1 2  3   4   5    6   7   8    9
    param = {8,10,500,3,[3,6],0.65,0,false,false};
end

dispfig = param{8};

%% Initialize vl_feat and obtain features
% addpath('vlfeat/toolbox');
addpath('../vlfeat/toolbox');
vl_setup

addpath(genpath('~/Projects/PeterKovesi/MatlabFns/'));

disp(['Reading target ' page1file ' .']);
page1 = imread(page1file);

disp(['Reading template ' page2file ' .']);
page2 = imread(page2file);

disp('Resizing template.');
page2 = imresize(page2,[size(page1,1) size(page1,2)]);

disp(['Registering ' page1file ' -> ' page2file]);


disp('Calculating SIFT features for target');
[fa, da] = vl_sift(single(rgb2gray(page1)) , 'PeakThresh', 0.01,'Magnif',param{1});
disp(['  Page 1 total features: ' num2str(length(fa))]);

disp('Calculating SIFT features for template');
[fb, db] = vl_sift(single(rgb2gray(page2)) , 'PeakThresh', 0.01,'Magnif',param{1});
disp(['  Page 2 total features: ' num2str(length(fb))]);

%% Feature matching
disp('Matching features');

matches = [];
scores = [];
prev_matches = [];
prev_scores = [];
limit = param{2};

% Start with strict matching and relax constraints until we roughly obtain
% a target number of matches. This solves the problem of some pages
% yielding a number of points while others yield very few at any one matching
% constraint.
while(length(matches) < param{3} && limit >= 2)
    prev_matches = matches;
    prev_scores = scores;
    [matches, scores] = vl_ubcmatch(da, db,limit);
    disp(['  ' num2str(limit) '. Matches: ' num2str(length(matches))]);
    limit = limit - 1;
end

% Overshot target; settle for fewer points
if( length(matches) >  param{3}+200 && limit < 10 )
    matches = prev_matches;
    scores = prev_scores;
end

if(dispfig), plotPages(page1,page2,fa,fb,matches,'Initial Matches'); end


% save temptest
% load temptest

%% Match filtering

disp('Filtering duplicate matches');
sel = filterMultipleMatches(matches);
[matches,scores,page1feat,page2feat] = reduceMatches(sel,matches,scores,fa,fb);

if(dispfig), plotPages(page1,page2,fa,fb,matches,'Duplicates Removed'); end

% disp('Filtering by homography');
% sel = filterHomography(matches,fa,fb,0.0001);
% [matches,scores,page1feat,page2feat] = reduceMatches(sel,matches,scores,fa,fb);
% 
% if(dispfig), plotPages(page1,page2,fa,fb,matches,'Homography'); end

disp('Filtering by fundamental matrix');
sel = filterFundamentalMatrix(matches,fa,fb,0.000001);
[matches,scores,page1feat,page2feat] = reduceMatches(sel,matches,scores,fa,fb);

if(dispfig), plotPages(page1,page2,fa,fb,matches,'Fundamental Matrix'); end

disp('Regressing to remove bad matches');
sel = regressMatches(page1feat,page2feat,matches,param{4});
[matches,scores,page1feat,page2feat] = reduceMatches(sel,matches,scores,fa,fb);

if(dispfig), plotPages(page1,page2,fa,fb,matches,'Regression'); end

disp('Neighborhood consistency check');
sel = localMatchFilter(page1feat,page2feat,matches,scores,size(page1),param{5},param{6},param{6});
[matches,scores,page1feat,page2feat] = reduceMatches(sel,matches,scores,fa,fb);

if(dispfig), plotPages(page1,page2,fa,fb,matches,'Neighborhood Consistency'); end

% Randomly sample a fixed number of matches
if param{7} > 0
    num = param{7};
    if size(matches,2) > num
        [matches,scores,page1feat,page2feat] = randomSel(matches,scores,fa,fb,num);
    end
    if(dispfig), plotPages(page1,page2,fa,fb,matches,'Random Sampling'); end
end

disp(['Final number of matches: ' num2str(size(matches,2)) ]);

% Allow the user to manually adjust the points
if param{9}
    [page2feat,page1feat] = cpselect(page2,page1,page2feat,page1feat,'Wait',true);
end

%% Registration

return;

if size(matches,2) > 2
    registered = registerImages(page1,page2,page1feat,page2feat);
    if(dispfig), figure; imshow(registered); title('Final Registration'); end
else
    disp('Not enough matches found');
    registered = [];
end


end


function plotPages(page1,page2,fa,fb,matches,title_str)
% PLOTPAGES     Utility function to display matches between pages.
xa = fa(1,matches(1,:));
% xb = fb(1,matches(2,:));
xb2 = fb(1,matches(2,:)) + size(page1,2) ;
ya = fa(2,matches(1,:));
yb = fb(2,matches(2,:));
figure;
imagesc(cat(2, page1(1:min(size(page1,1),size(page2,1)),1:min(size(page1,2),size(page2,2)),:), page2(1:min(size(page1,1),size(page2,1)),1:min(size(page1,2),size(page2,2)),:)));
hold on ;
h = line([xa ; xb2], [ya ; yb]);
set(h,'linewidth', 2, 'color', 'b');
vl_plotframe(fa(:,matches(1,:)));
fb(1,:) = fb(1,:) + size(page1,2);
vl_plotframe(fb(:,matches(2,:)));
axis equal;
axis off;
if(nargin > 5)
    title(title_str);
end
end

function [matches,scores,page1feat,page2feat] = reduceMatches(sel,matches,scores,fa,fb)
% REDUCEMATCHES     Utility function to update all relevant variables with
%                   a new selection.
matches = matches(:,sel);
scores = scores(sel);
xa = fa(1,matches(1,:));
xb = fb(1,matches(2,:));
ya = fa(2,matches(1,:));
yb = fb(2,matches(2,:));
page1feat = [xa' ya'];
page2feat = [xb' yb'];
end

function [matches,scores,page1feat,page2feat] = randomSel(matches,scores,fa,fb,num)
%RANDOMSEL     Randomly sample remaining matches.

    % randn('state',0) ;
    % rand('state',0) ;
perm = randperm(size(matches,2));
sel = perm(1:num);
matches = matches(:,sel);
scores = scores(sel);
xa = fa(1,matches(1,:));
xb = fb(1,matches(2,:));
ya = fa(2,matches(1,:));
yb = fb(2,matches(2,:));
page1feat = [xa' ya'];
page2feat = [xb' yb'];
end


function final_matches = filterMultipleMatches(matches)
% FILTERMULTIPLEMATCHES     Remove matches that have features in common
counts = zeros(max(matches(2,:)),1);

% Create bins that store the "use counts" for each feature
for i = 1:size(matches,2)
    counts(matches(2,i)) = counts(matches(2,i)) + 1;
end

final_matches = [];

% Only add back matches that have a "use count" of 1
for i = 1:size(matches,2)
    if counts(matches(2,i)) == 1
        final_matches = [final_matches  i]; % matches(:,i) ];
    end
end

disp(['Removed ' num2str(size(matches,2)-size(final_matches,2)) ' duplicates.' ]);

end

function final_matches = filterHomography(matches,fa,fb,t)

m1 = [fa(2,matches(1,:)); fa(1,matches(1,:))];
m2 = [fb(2,matches(2,:)); fb(1,matches(2,:))];

% Display putative matches
% figure; imshow(page1); hold on;
% for n = 1:length(m1);
%     line([m1(2,n) m2(2,n)], [m1(1,n) m2(1,n)])
% end

x1 = [m1(2,:); m1(1,:); ones(1,length(m1))];
x2 = [m2(2,:); m2(1,:); ones(1,length(m1))];

% t = 0.0001;  % Distance threshold for deciding outliers
[H, inliers] = ransacfithomography(x1, x2, t);

fprintf('Number of inliers was %d (%d%%) \n', ...
    length(inliers),round(100*length(inliers)/length(m1)))
fprintf('Number of putative matches was %d \n', length(m1))

% figure; imshow(page1); hold on;
% plot(m1(2,inliers),m1(1,inliers),'r+');
% plot(m2(2,inliers),m2(1,inliers),'g+');
% for n = inliers
%     line([m1(2,n) m2(2,n)], [m1(1,n) m2(1,n)],'color',[0 0 1])
% end

final_matches = inliers;

end

function final_matches = filterFundamentalMatrix(matches,fa,fb,t)

m1 = [fa(2,matches(1,:)); fa(1,matches(1,:))];
m2 = [fb(2,matches(2,:)); fb(1,matches(2,:))];

x1 = [m1(2,:); m1(1,:); ones(1,length(m1))];
x2 = [m2(2,:); m2(1,:); ones(1,length(m1))];

% t = 0.0001;  % Distance threshold for deciding outliers
[F, inliers] = ransacfitfundmatrix(x1, x2, t, 1);

fprintf('Number of inliers was %d (%d%%) \n', ...
    length(inliers),round(100*length(inliers)/length(m1)))
fprintf('Number of putative matches was %d \n', length(m1))

final_matches = inliers;

end