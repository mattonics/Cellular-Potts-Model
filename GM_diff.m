function d = GM_diff(activities, grid, u_I, v_I)
    [x,y] = ind2sub(size(activities), u_I);
    [nx, ny] = meshgrid(x-1:x+1, y-1:y+1);
    n = [nx(:) ny(:)]; 
    n(n(:,1) == 0 | n(:,1) == size(activities,1)+1 | n(:,2) == 0 | n(:,2) == size(activities,2)+1,:) = [];
    inds = sub2ind(size(activities),n(:,1),n(:,2));
    g1 = geomean(activities(inds(grid(inds) == grid(u_I))));
    
    [x,y] = ind2sub(size(activities), v_I);
    [nx, ny] = meshgrid(x-1:x+1, y-1:y+1);
    n = [nx(:) ny(:)]; 
    n(n(:,1) == 0 | n(:,1) == size(activities,1)+1 | n(:,2) == 0 | n(:,2) == size(activities,2)+1,:) = [];
    inds = sub2ind(size(activities),n(:,1),n(:,2));
    g2 = geomean(activities(inds(grid(inds) == grid(v_I))));
    
    d = g1 - g2;
end