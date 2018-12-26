u = VideoWriter('C:\Users\Matt\Desktop\CPM_singlecell.avi','Motion JPEG AVI');
open(u);

Max_Act = 80;
lambda_Act = 200;

grid = zeros(100,100);
grid(30:60,40:70) = 1;
J_Cell_Medium = 20;
J_Cell_Cell = 40;
volume_constraint = 50;
SA_constraint = 2;
Temperature = 20;

A0 = 200;
L0 = 180;

p = perim(grid,1);
volume = sum(grid(:));
activities = zeros(size(grid));
gcf = figure(1);

for j = 1:100
    j
    randomCells = randi([1 numel(grid)],[1 round(numel(grid))]);
    [x,y] = ind2sub(size(grid),randomCells); 

    for i = 1:length(randomCells)
        [nx, ny] = meshgrid(x(i)-1:x(i)+1, y(i)-1:y(i)+1);
        n = [nx(:) ny(:)]; 
        n(n(:,1) == 0 | n(:,1) == size(grid,1)+1 | n(:,2) == 0 | n(:,2) == size(grid,2)+1,:) = [];
        inds = sub2ind(size(grid),n(:,1),n(:,2));
        all_neighbors = grid(inds);
        ri = randi([1 length(all_neighbors)]);
        rand_neighbor = all_neighbors(ri);
        ii = inds(ri);
        if ( grid(x(i), y(i) ) ~= rand_neighbor )
            H = (p - L0)^2*SA_constraint+(volume-A0)^2*volume_constraint;        
            newgrid = grid;
            newgrid(x(i), y(i) ) = rand_neighbor;
            newp = perim(newgrid,1);
            if ( rand_neighbor == 0 )
                newv = volume - 1;
            else
                newv = volume + 1;
            end
            oldH_J = H_J(grid, J_Cell_Cell, J_Cell_Medium, x(i), y(i));
            newH_J = H_J(newgrid, J_Cell_Cell, J_Cell_Medium, x(i), y(i));
            newH = (newp - L0)^2*SA_constraint+(newv-A0)^2*volume_constraint;
            deltaH2 = newH - H - newH_J + oldH_J;
            
            if ( test_connectivity(grid, sub2ind(size(grid),x(i),y(i))) > 2)
                deltaH2 = deltaH2 + 3000;
            end
            pr2 = exp(-deltaH2/Temperature);
            if ( pr2 > rand() )
                grid = newgrid;
                p = newp;
                volume = newv;
                H = newH;
            end            
        end
    end
    activities = zeros(size(grid));
    activities(logical(grid)) = Max_Act;    
    subplot(1,2,1)
    imshow(grid);
    colormap('hot');
    subplot(1,2,2)
    imshow(activities,[0 Max_Act]);
    colormap('hot');
    colorbar
    truesize(gcf,[200 200])

    writeVideo(u,getframe(gcf));
end
imshow(grid);

for j = 1:300
    j
    randomCells = randi([1 numel(grid)],[1 round(numel(grid))]);
    [x,y] = ind2sub(size(grid),randomCells); 

    for i = 1:length(randomCells)
        [nx, ny] = meshgrid(x(i)-1:x(i)+1, y(i)-1:y(i)+1);
        n = [nx(:) ny(:)]; 
        n(n(:,1) == 0 | n(:,1) == size(grid,1)+1 | n(:,2) == 0 | n(:,2) == size(grid,2)+1,:) = [];
        inds = sub2ind(size(grid),n(:,1),n(:,2));
        all_neighbors = grid(inds);
        ri = randi([1 length(all_neighbors)]);
        rand_neighbor = all_neighbors(ri);
        ii = inds(ri);
        if ( grid(x(i), y(i) ) ~= rand_neighbor )
            H = (p - L0)^2*SA_constraint+(volume-A0)^2*volume_constraint;        
            newgrid = grid;
            newgrid(x(i), y(i) ) = rand_neighbor;
            newp = perim(newgrid,1);
            if ( rand_neighbor == 0 )
                newv = volume - 1;
            else
                newv = volume + 1;
            end
            oldH_J = H_J(grid, J_Cell_Cell, J_Cell_Medium, x(i), y(i));
            newH_J = H_J(newgrid, J_Cell_Cell, J_Cell_Medium, x(i), y(i));
            newH = (newp - L0)^2*SA_constraint+(newv-A0)^2*volume_constraint;
            deltaH2 = newH - H - newH_J + oldH_J;
            deltaH2 = deltaH2 - lambda_Act/Max_Act * GM_diff(activities, grid, ii, sub2ind(size(grid),x(i),y(i)));
            
            if ( test_connectivity(grid, sub2ind(size(grid),x(i),y(i))) > 2)
                deltaH2 = deltaH2 + 3000;
            end
            pr2 = exp(-deltaH2/Temperature);
            if ( pr2 > rand() )
                if ( rand_neighbor == 0 )
                    activities(x(i), y(i)) = 0;
                elseif ( rand_neighbor == 1 )
                    activities(x(i),y(i)) = Max_Act;
                end
                grid = newgrid;
                p = newp;
                volume = newv;
                H = newH;
            end            
        end
    end
    subplot(1,2,1)
    imshow(grid);
    colormap('hot');
    subplot(1,2,2)
    imshow(activities,[0 Max_Act]);
    colormap('hot');
    colorbar
    truesize(gcf,[200 200])

    writeVideo(u,getframe(gcf));
    activities(activities ~= 0) = activities(activities ~= 0) - 1;
end
close(u);