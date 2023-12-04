function volume = sphereIntersection(circleCenters,radius)
    % Used to estimate the SOZ, resected, and resected SOZ volume by
    % inflating each channel by a sphere and computing the total volume of
    % the resulting mask

    if size(circleCenters,1) == 0
        volume = 0;
        return
    end
    % Inflating contacts by radius
    coords_new = [];

    N = 100;
    thetavec = linspace(0,pi,N);
    phivec = linspace(0,2*pi,2*N);
    [th, ph] = meshgrid(thetavec,phivec);
    R = radius*ones(size(th)); % should be your R(theta,phi) surface in general
    
    x = R.*sin(th).*cos(ph); x = reshape(x,[numel(x),1]);
    y = R.*sin(th).*sin(ph); y = reshape(y,[numel(y),1]);
    z = R.*cos(th); z = reshape(z,[numel(z),1]);

    coords_new = [];
    for i = 1:size(circleCenters,1)
        coords = circleCenters(i,:);
        
        % Inflate by radius
        coords = [x + coords(1) y + coords(2) z + coords(3)]; 
        % tmp = 

        coords_new = [coords_new; coords];
    end

    % Find the minimum and maximum x and y coordinates for the circles
    minX = min(coords_new(:, 1) - radius);
    maxX = max(coords_new(:, 1) + radius);
    minY = min(coords_new(:, 2) - radius);
    maxY = max(coords_new(:, 2) + radius);
    minZ = min(coords_new(:, 3) - radius);
    maxZ = max(coords_new(:, 3) + radius);
    
    % Calculate the width and height of the square
    width = maxX - minX;
    height = maxY - minY;
    depth = maxZ - minZ;
    
    
    % Generate a grid of points within the volume
    x = (minX:0.5:maxX)';
    y = (minY:0.5:maxY)';
    z = (minZ:0.5:maxZ)';
    
    % Create a grid of points using meshgrid
    [xx, yy, zz] = meshgrid(x, y, z);
    
    b = 0&1;
    for i = 1:size(circleCenters,1)
        c=(xx-circleCenters(i,1)).^2+(yy-circleCenters(i,2)).^2+(zz-circleCenters(i,3)).^2; 
        b = b | ( c <= radius^2);
    end
 
    
    volumeSquare = width*height*depth;
    volume = volumeSquare*sum(b,'all')/numel(b);
end