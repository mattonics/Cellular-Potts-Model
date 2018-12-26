videoWrite = 1;

if (videoWrite)
    v = VideoWriter('C:\Users\Matt\Desktop\CPM_manycells.avi','Motion JPEG AVI');
    open(v);
end

Max_Act = 80;
lambda_Act = 4000;

p = 0.3;
grid = zeros(200,200);
i = 1;
cell_size = 20;
for j = 1:size(grid,1)/cell_size
    for k = 1:size(grid,2)/cell_size
        if ( rand() < p )
           grid((j-1)*cell_size+1:j*cell_size,(k-1)*cell_size+1:k*cell_size) = i;
           i = i+1;
        end
    end
end

J_Cell_Medium = 20;
J_Cell_Cell = 40;
volume_constraint = 50;
SA_constraint = 2;
Temperature = 20;


A0 = 200;
L0 = 180;

p = perim_mex(grid,[]);
p = p(2:end);
newp = p;
volume = vol(grid);
newv = volume;
gca = figure(1);
for j = 1:100
    j
    randomCells = randi([1 numel(grid)],[1 numel(grid)]);
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
            allIs = [grid(x(i), y(i) ) rand_neighbor];
            allIs(allIs == 0) = [];
            Hp = sum((p(allIs)-L0).^2*SA_constraint);
            Hv = sum((volume(allIs)-A0).^2*volume_constraint);
            H = Hp+Hv;
            newgrid = grid;
            newgrid(x(i), y(i) ) = rand_neighbor;
            prev = grid(x(i), y(i) );
            
            np = perim_mex(newgrid,[prev rand_neighbor]);
            newp = p;
            if (np(1) ~= -1)
                newp(prev) = np(1);
            end
            if (np(2) ~= -1)
                newp(rand_neighbor) = np(2);
            end
            newv = volume;
            if (prev == 0)
                newv(rand_neighbor) = newv(rand_neighbor)+1;
            elseif (rand_neighbor == 0)
                newv(prev) = newv(prev) - 1;
            else
                newv(prev) = newv(prev) - 1;
                newv(rand_neighbor) = newv(rand_neighbor)+1;
            end

            newH_J_mex = H_J_mex(newgrid, J_Cell_Cell, J_Cell_Medium, x(i), y(i));
            oldH_J_mex = H_J_mex(grid, J_Cell_Cell, J_Cell_Medium, x(i), y(i));
            newHp = sum(SA_constraint*(newp(allIs)-L0).^2);
            newHv = sum(volume_constraint*(newv(allIs)-A0).^2);
            newH = newHp+newHv;
            deltaH2 = newH - H - newH_J_mex + oldH_J_mex;
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
    if (videoWrite)
        activities = zeros(size(grid));
        activities(logical(grid)) = Max_Act;    
        s1 = subplot(1,2,1);
        imshow(grid/max(grid(:)));
        colormap(s1, 'colorcube');
        s2 = subplot(1,2,2);
        imshow(activities,[0 Max_Act]);
        colormap(s2, 'hot');
        colorbar
        truesize(gca,[200 200]);
        writeVideo(v,getframe(gca));
    end
end

activities(grid ~= 0) = Max_Act;

for j = 1:300
    j
    randomCells = randi([1 numel(grid)],[1 numel(grid)]);
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
            allIs = [grid(x(i), y(i) ) rand_neighbor];
            allIs(allIs == 0) = [];
            Hp = sum((p(allIs)-L0).^2*SA_constraint);
            Hv = sum((volume(allIs)-A0).^2*volume_constraint);
            H = Hp+Hv;
            newgrid = grid;
            newgrid(x(i), y(i) ) = rand_neighbor;
            prev = grid(x(i), y(i) );

            newp = p;
            np = perim_mex(newgrid,[prev rand_neighbor]);
            if (np(1) ~= -1)
                newp(prev) = np(1);
            end
            if (np(2) ~= -1)
                newp(rand_neighbor) = np(2);
            end

            newv = volume;
            if (prev == 0)
                newv(rand_neighbor) = newv(rand_neighbor)+1;
            elseif (rand_neighbor == 0)
                newv(prev) = newv(prev) - 1;
            else
                newv(prev) = newv(prev) - 1;
                newv(rand_neighbor) = newv(rand_neighbor)+1;
            end
            newH_J_mex = H_J_mex(newgrid, J_Cell_Cell, J_Cell_Medium, x(i), y(i));
            oldH_J_mex = H_J_mex(grid, J_Cell_Cell, J_Cell_Medium, x(i), y(i));
            newHp = sum(SA_constraint*(newp(allIs)-L0).^2);
            newHv = sum(volume_constraint*(newv(allIs)-A0).^2);
            newH = newHp+newHv;
            deltaH2 = newH - H - newH_J_mex + oldH_J_mex;
            deltaH2 = deltaH2 - lambda_Act/Max_Act * GM_diff(activities, grid, ii, sub2ind(size(grid),x(i),y(i)));
            if ( test_connectivity(grid, sub2ind(size(grid),x(i),y(i))) > 2)
                deltaH2 = deltaH2 + 3000;
            end
            pr2 = exp(-deltaH2/Temperature);
            if ( pr2 > rand() )
                if ( rand_neighbor == 0 )
                    activities(x(i), y(i)) = 0;
                elseif ( rand_neighbor ~= 0 )
                    activities(x(i),y(i)) = Max_Act;
                end
                grid = newgrid;
                p = newp;
                volume = newv;
                H = newH;
            end
        end
    end
    
    if (videoWrite)
        s1 = subplot(1,2,1);
        imshow(grid/max(grid(:)));
        colormap(s1, 'colorcube');
        s2 = subplot(1,2,2);
        imshow(activities,[0 Max_Act]);
        colormap(s2, 'hot');
        colorbar
        truesize(gca,[200 200]);
        writeVideo(v,getframe(gca));
    end
    
    activities(activities ~= 0) = activities(activities ~= 0) - 1;

end
if (videoWrite)
    close(v);
end
imshow(grid);