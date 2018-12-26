function sum = test_connectivity(grid, i)
    [x,y] = ind2sub(size(grid),i);
    x_i = grid(i);

    [nx, ny] = meshgrid(x-1:x+1, y-1:y+1);
    n = [nx(:) ny(:)]; 
    clockwise = [n(1:3,:);n(6,:);n(9,:);n(8,:);n(7,:);n(4,:)];
    values = zeros(8,1);
    temp = clockwise;
    values(clockwise(:,1) == 0 | clockwise(:,1) == size(grid,1)+1 | clockwise(:,2) == 0 | clockwise(:,2) == size(grid,2)+1) = 0;
    temp(clockwise(:,1) == 0 | clockwise(:,1) == size(grid,1)+1 | clockwise(:,2) == 0 | clockwise(:,2) == size(grid,2)+1,:) = [];
    goodvalues = grid(sub2ind(size(grid),temp(:,1),temp(:,2)));
    values(~(clockwise(:,1) == 0 | clockwise(:,1) == size(grid,1)+1 | clockwise(:,2) == 0 | clockwise(:,2) == size(grid,2)+1)) = goodvalues;
    sum = 0;
    for a = 1:length(values)
        j = a - 1;
        if (j == 0)
            j = length(values);
        end
        k = a+1;
        if (k == length(values)+1)
            k = 1;
        end
        sum = sum + (values(a) == x_i)*(2 - (values(k) == x_i) - (values(j) == x_i));
    end
end