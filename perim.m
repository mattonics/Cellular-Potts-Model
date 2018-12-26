function u = perim(grid, updates)

if isempty(updates)
    updates = unique(grid(:));
end

u = updates;
updates = updates(updates > 0);

s = zeros(length(updates),1);
for j = 1:length(updates)
    g = grid;
    g(g ~= updates(j)) = 0;
    g(g == updates(j)) = 1;
    [x,y] = ind2sub(size(g),find(g));
    s(j) = 0;
    for i = 1:length(x)
        [nx,ny] = meshgrid(x(i)-1:x(i)+1,y(i)-1:y(i)+1);
        n = [nx(:) ny(:)];
        n(n(:,1) == 0 | n(:,1) == size(g,1)+1 | n(:,2) == 0 | n(:,2) == size(grid,2)+1,:) = [];
        n_I = sub2ind(size(grid),n(:,1),n(:,2));
        n_I(grid(n_I) < 0) = [];
        all_neighbors = g(n_I);
        if ( g(x(i),y(i)) == 1 )
            s(j) = s(j) + sum(~all_neighbors);
        end
    end
end

u(u > 0) = s;
u(u == 0) = -1;

end