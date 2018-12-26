function s = H_J(grid, J_Cell_Cell, J_Cell_Medium, x, y)

[nx, ny] = meshgrid(x-1:x+1, y-1:y+1);
n = [nx(:) ny(:)]; 
n(n(:,1) <= 0 | n(:,1) >= size(grid,1)+1 | n(:,2) <= 0 | n(:,2) >= size(grid,2)+1,:) = [];
n(grid(sub2ind(size(grid),n(:,1),n(:,2)) < 0),:) = [];

x = n(:,1); y = n(:,2);
s = 0;
for i = 1:length(x)
    [nx, ny] = meshgrid(x(i)-1:x(i)+1, y(i)-1:y(i)+1);
    n = [nx(:) ny(:)]; 
    n(n(:,1) == 0 | n(:,1) == size(grid,1)+1 | n(:,2) == 0 | n(:,2) == size(grid,2)+1,:) = [];
    inds = sub2ind(size(grid),n(:,1),n(:,2));
    all_neighbors = grid(inds);
    if ( grid(x(i),y(i)) == 0 )
        s = s + sum((all_neighbors > 0))*J_Cell_Medium;
    elseif ( grid(x(i),y(i)) > 0 )
        s = s + sum(all_neighbors == 0)*J_Cell_Medium;
        s = s + sum(all_neighbors > 0 & all_neighbors ~= grid(x(i),y(i)))*J_Cell_Cell;
    end
end
end