function final_matches = localMatchFilter(page1feat,page2feat, matches, scores, dims, slices, xthresh, ythresh)
% LOCALMATCHFILTER     Score each match by doing a statistical jackknife on
% each match and its local neighborhood (of 20 nearest neighbors). A
% percentage of the worst scores within each grid cell of the image are
% removed from each cell.

%% Initialization

NNNum = 30;
exclNum = 11;

verbose = false;

n=length(matches);

vx = page1feat(:,1);
vy = page1feat(:,2);
ux = page2feat(:,1);
uy = page2feat(:,2);

scores = scores - min(scores);
scores = scores ./ max(scores);

border_pad = [ceil(dims(2)/15) ceil(dims(1)/15)];


differences = [];
vars = [];

%% Compute score for each match and store in "differences"
for i = 1:n
    % Special cases to retain order of matches found by knnsearch
    if(i == 1)
        idx = knnsearch([vx(i) vy(i)],page1feat(2:end,:),NNNum);
    elseif(i == n)
        idx = knnsearch([vx(i) vy(i)],page1feat(1:end-1,:),NNNum);
    else
        idx = knnsearch([vx(i) vy(i)],[page1feat(1:i-1,:) ; page1feat(i+1:end,:)],NNNum);
    end
    
    % Remove fixed number of "similar" matches to the match under
    % consideration in case there are clusters of bad matches that
    % happen to map under the same general transform
    [temp_sort,exclude_idx] = sort(sqrt(((page1feat(idx,1)-page2feat(idx,1))-(page1feat(i,1)-page2feat(i,1))).^2 + ((page1feat(idx,2)-page2feat(idx,2))-(page1feat(i,2)-page2feat(i,2))).^2 ),'descend');
    idx(exclude_idx(1:exclNum)) = [];
    
    % Compute average change in distance /without/ center point
    page1pts = [page1feat(idx,1) page1feat(idx,2)];
    page2pts = [page2feat(idx,1) page2feat(idx,2)];
    dist = abs(page1pts-page2pts);    
    avgdist = mean(dist,1);
    
    % Compute average change in distance /with/ center point
    page1pts2 = [page1feat(idx,1) page1feat(idx,2) ; page1feat(i,1) page1feat(i,2)];
    page2pts2 = [page2feat(idx,1) page2feat(idx,2) ; page2feat(i,1) page2feat(i,2)];
    dist2 = abs(page1pts2-page2pts2);
    avgdist2 = mean(dist2,1);
    
    % Score is given by the difference of average distances (i.e. change
    % between resampled averages).
    differences = [differences; abs(avgdist-avgdist2)*scores(i)];
end

%% Divide page into a grid and remove bad matches from each cell

final_matches = [];
reconsider_matches = [];

% slices = [3,6];

% Determine the size of each chunk
xchunks = floor(dims(2)/slices(1));
ychunks = floor(dims(1)/slices(2));

% Temp variable to keep track of the previous chunk size
prev = [1,1];

% X chunk
for k = 1:slices(1)    
    % Y chunk
    for l = 1:slices(2)
        
        if(verbose), disp(['k=' num2str(k) ', l=' num2str(l)]); end
        
        % Handle border cases
        if k == 1
            xlims = 1:xchunks;
        elseif k == slices(1)
            xlims = prev(1):dims(2);
        else 
            xlims = xchunks*(k-1):xchunks*k;
            prev(1) = xchunks*k;
        end
        
        if l == 1
            ylims = 1:ychunks;
        elseif l == slices(2)
            ylims = prev(2):dims(1);
        else 
            ylims = ychunks*(l-1):ychunks*l;
            prev(2) = ychunks*l;
        end
        
        if(verbose), disp(['From [' num2str(xlims(1)) ',' num2str(xlims(length(xlims))) '] to [' num2str(ylims(1)) ',' num2str(ylims(length(ylims))) ']' ]); end
        
        % Get the vectors of xs and ys that are within the box,
        new_idx = find(ismember(round(ux),xlims) & ismember(round(uy),ylims)==1);
        
        if length(new_idx) < 2
            disp('Too few matches in region.');
            if length(new_idx) == 1
                % We have to trust that the regression removed lonely outliers
                reconsider_matches = [reconsider_matches new_idx];
            end
            continue;
        end
        
        if(verbose), disp(['In square: ' num2str(length(new_idx)) ]); end
        
        % Sort and select outlier point based on the given threshold
        sorted_differencesx = sort(differences(new_idx,1));
        sorted_differencesy = sort(differences(new_idx,2));
        
        outlier_val = [];
        
        outlier_val(1) = sorted_differencesx(min(max(round(xthresh*length(sorted_differencesx)),1),length(sorted_differencesx)-1));
        outlier_val(2) = sorted_differencesy(min(max(round(ythresh*length(sorted_differencesy)),1),length(sorted_differencesy)-1));
        
        if(verbose), new_count = 0; end

        % Iterate through all matches in cell that are below the threshold        
        for i = new_idx'
            if(differences(i,1) <= outlier_val(1) && differences(i,2) <= outlier_val(2))
                % Enforce an additional condition that the matched points
                % not be too close to the border.
                if ( ~( ((uy(i) > dims(1)-border_pad(2)) || (uy(i) < border_pad(2))) && ((vy(i) > dims(1)-border_pad(2)) || (vy(i) < border_pad(2))) ) && ...
                    ~( ((ux(i) > dims(2)-border_pad(1)) || (ux(i) < border_pad(1))) || ((vx(i) > dims(2)-border_pad(1)) || (vx(i) < border_pad(1))) ) )
                    final_matches = [final_matches i];
                end
                if(verbose), new_count = new_count + 1; end
            end
        end
        
        if(verbose), disp(['Removed ' num2str(length(sorted_differencesx) - new_count)]); end
        
    end
end



disp(['Matches removed: ' num2str(length(matches)-length(final_matches))]);

end
