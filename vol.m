function s = vol(grid)

types = unique(grid);
types = types(types > 0);

for j = 1:length(types)
    g = grid;
    g(g ~= types(j)) = 0;
    g(g == types(j)) = 1;
    s(j) = sum(g(:));
end

end