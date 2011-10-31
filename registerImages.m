function registered = registerImages(base_file,unregistered_file,base_landmarks,unregistered_landmarks)

if nargin < 3

fid1 = fopen([unregistered_file '_landmarks']);
if(fid1 ~= -1)
    fclose(fid1);
    unregistered_landmarks = load([unregistered_file '_landmarks']);
else
    unregistered_landmarks = [];
end

fid2 = fopen([base_file '_landmarks']);
if(fid2 ~= -1)
    fclose(fid2);
    base_landmarks = load([base_file '_landmarks']);
else
    base_landmarks = [];
end

%read in files
unregistered = imread(unregistered_file);
base = imread(base_file);

%provide tool for user to select matching landmarks
if(isempty(base_landmarks)==0 || isempty(unregistered_landmarks) == 0)
    [unregistered_landmarks, base_landmarks] = cpselect(unregistered,base,unregistered_landmarks,base_landmarks,'Wait',true);
else
    [unregistered_landmarks, base_landmarks] = cpselect(unregistered,base,'Wait',true);
end
dlmwrite([unregistered_file '_landmarks'],unregistered_landmarks,' ');
dlmwrite([base_file '_landmarks'],base_landmarks,' ');

else
    unregistered = unregistered_file;
    base = base_file;
end

% figure; imshow(unregistered); hold on; plot(unregistered_landmarks(:,1),unregistered_landmarks(:,2),'r.'); 
% print('-dtiff',[unregistered_file '_POINTS']); close all;
% figure; imshow(base); hold on; plot(base_landmarks(:,1),base_landmarks(:,2),'g.'); 
% print('-dtiff',[base_file '_POINTS']); close all;


%calculate transorm
% tform = findTransform(unregistered_landmarks,base_landmarks);
base_landmarks = fliplr(base_landmarks);
unregistered_landmarks=fliplr(unregistered_landmarks);
tform = findTransform(base_landmarks,unregistered_landmarks);

%transform image
registered = base;
% registered(:,:,:) = 0;
[registered,tlandmarks] = tpsTransformImage(base,unregistered,tform,base_landmarks);
fliplr(tlandmarks);

% figure; imshow(registered);
% hold on; h = imshow(base,gray(256));
% set(h,'AlphaData',0.6);
% % print('-dpng',[unregistered_file '_OUT']);
% close all;

%write result

%imwrite(registered,[unregistered_file '_OUT'],'tiff');

% figure; imshow(registered); 

% hold on; plot(tlandmarks(:,1),tlandmarks(:,2),'b.'); 
% print('-dtiff',[unregistered_file '_RES']); close all;

end
