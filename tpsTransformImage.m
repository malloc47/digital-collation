function [registered,tlandmarks] = tpsTransformImage(base,unregistered,tform,landmarks)
% TPSTRANSFORMIMAGE     Compute the transform of every point in the image
%                       with the given transform (tform) resulting in a 
%                       registered image as output.
m = length(tform);
n = length(landmarks);
[h,w,c] = size(base);
registered = base;%zeros(h,w,c);

[X Y C] = size(unregistered);

tpsW = tform(1:n,:);
tpsT = tform(n+1,:);
tpsA = tform(n+2:n+3,:);

tlandmarks = zeros(size(landmarks));

for i=1:n
    val = tpsT;
    val = val+landmarks(i,:)*tpsA(:,:);
    for k=1:n
        r = (landmarks(i,:) - landmarks(k,:)) * (landmarks(i,:) - landmarks(k,:))';
        r = sqrt(r);
        if(r > 0.0)

            val(1) = val(1) + tpsW(k,1)*r*r*log(r);
            val(2) = val(2) + tpsW(k,2)*r*r*log(r);
        end
    end
    tlandmarks(i,:) = val;
end

for i=1:h
    for j = 1:w
        
        val = tpsT;
        val = val+[i j]*tpsA(:,:);
        
        for k=1:n
            r = ([i j] - landmarks(k,:)) * ([i j] - landmarks(k,:))';
            r = sqrt(r);
            if(r > 0.0)
                val(1) = val(1) + tpsW(k,1)*r*r*log(r);
                val(2) = val(2) + tpsW(k,2)*r*r*log(r);
            end
        end
        x = round(val(1));
        y = round(val(2));

        if (x > X || y > Y || x < 1 || y < 1)
            registered(i,j,:) = 0;
        else
            registered(i,j,:) = unregistered(x,y,:);
        end

    end
end

end