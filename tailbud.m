tb = zeros(474,402);

borderValue = -1;

tb(292:end,92) = borderValue;
tb(292:end,111) = borderValue;
tb(292:end,292) = borderValue;
tb(292:end,311) = borderValue;
tb(292,92:111) = borderValue;
tb(292, 292:311) = borderValue;
tb(292:end,92:111) = borderValue;
tb(292:end, 292:311) = borderValue;

center = [202 201.5];
th = linspace( pi, 2*pi, 1000);
R = 200.5;  %or whatever radius you want
x = R*cos(th) + center(2);
y = R*sin(th) + center(1);

tb(sub2ind(size(tb),round(y)',round(x)')) = borderValue;
tb(sub2ind(size(tb),round(y)'-1,round(x)')) = borderValue;

tb(end,:) = borderValue;

imshow(tb)
for i = 1:size(tb,2)
    col = tb(:,i);
    j = find(col,1);
    if (j ~= 1)
        tb(1:(j-1),i) = borderValue;
    end
end

p = 0.3;
grid = zeros(size(tb));
i = 1;
cell_size = 13;
for j = 1:size(tb,1)/cell_size
    for k = 1:size(tb,2)/cell_size
        if ( rand() < p )
           grid((j-1)*cell_size+1:j*cell_size,(k-1)*cell_size+1:k*cell_size) = i;
           i = i+1;
        end
    end
end

currcell = max(grid(:)) + 1;
for j = 1:15
    for k = 1:14
        grid(end-12*j+1:end-12*(j-1),111+12*(k-1):111+12*k - 1) = currcell;
        currcell = currcell + 1;
    end
end


grid(tb == borderValue) = -1;
imshow(grid);

cells = unique(grid(:));
cells(cells < 1) = [];
for i = 1:length(cells)
    grid(grid == cells(i)) = i;
end

J_Cell_Medium = 20;
J_Cell_Cell = 40;
volume_constraint = 50;
SA_constraint = 2;
Temperature = 20;

videoWrite = 1;

if (videoWrite)
    v = VideoWriter('C:\Users\Matt\Desktop\CPM_tailbud.avi','Motion JPEG AVI');
    open(v);
end

Max_Act = 60;
lambda_Act = 4000;

A0 = 180;
L0 = 162;

p = perim(grid,[]);
p = p(p ~= -1);
newp = p;
volume = vol(grid);
newv = volume;
% gca = figure(1);
% for j = 1:20
%     j
%     tailbudI = find(grid ~= borderValue);
%     randomCells = randi([1 numel(tailbudI)],[1 numel(tailbudI)]);
%     randomCells = tailbudI(randomCells);
%     [x,y] = ind2sub(size(grid),randomCells); 
% 
%     for i = 1:length(randomCells)
%         [nx, ny] = meshgrid(x(i)-1:x(i)+1, y(i)-1:y(i)+1);
%         n = [nx(:) ny(:)]; 
%         n(n(:,1) == 0 | n(:,1) == size(grid,1)+1 | n(:,2) == 0 | n(:,2) == size(grid,2)+1,:) = [];
%         inds = sub2ind(size(grid),n(:,1),n(:,2));
%         all_neighbors = grid(inds);
%         ri = randi([1 length(all_neighbors)]);
%         rand_neighbor = all_neighbors(ri);
%         ii = inds(ri);
%         if ( grid(x(i), y(i) ) ~= rand_neighbor && grid(x(i),y(i)) ~= -1 && rand_neighbor ~= -1)
%             allIs = [grid(x(i), y(i) ) rand_neighbor];
%             allIs(allIs == 0) = [];
%             Hp = sum((p(allIs)-L0).^2*SA_constraint);
%             Hv = sum((volume(allIs)-A0).^2*volume_constraint);
%             H = Hp+Hv;
%             newgrid = grid;
%             newgrid(x(i), y(i) ) = rand_neighbor;
%             prev = grid(x(i), y(i) );
%             np = perim(newgrid,[prev rand_neighbor]);
%             newp = p;
%             if (np(1) ~= -1)
%                 newp(prev) = np(1);
%             end
%             if (np(2) ~= -1)
%                 newp(rand_neighbor) = np(2);
%             end
%             newv = volume;
%             if (prev == 0)
%                 newv(rand_neighbor) = newv(rand_neighbor)+1;
%             elseif (rand_neighbor == 0)
%                 newv(prev) = newv(prev) - 1;
%             else
%                 newv(prev) = newv(prev) - 1;
%                 newv(rand_neighbor) = newv(rand_neighbor)+1;
%             end
% 
%             newH_J = H_J(newgrid, J_Cell_Cell, J_Cell_Medium, x(i), y(i));
%             oldH_J = H_J(grid, J_Cell_Cell, J_Cell_Medium, x(i), y(i));
%             newHp = sum(SA_constraint*(newp(allIs)-L0).^2);
%             newHv = sum(volume_constraint*(newv(allIs)-A0).^2);
%             newH = newHp+newHv;
%             deltaH2 = newH - H - newH_J + oldH_J;
%             if ( test_connectivity(grid, sub2ind(size(grid),x(i),y(i))) > 2)
%                 deltaH2 = deltaH2 + 3000;
%             end
%             pr2 = exp(-deltaH2/Temperature);
%             if ( pr2 > rand() )
%                 
%                 grid = newgrid;
%                 p = newp;
%                 volume = newv;
%                 H = newH;
%             end
%         end
%     end
%     if (videoWrite)
%         s1 = subplot(1,2,1);
%         g = grid;
%         g(g < 0) = 0;
%         activities = zeros(size(g));
%         activities(g > 0) = Max_Act;   
%         imshow(g/max(g(:)));
%         colormap(s1, 'colorcube');
%         s2 = subplot(1,2,2);
%         imshow(activities,[0 Max_Act]);
%         colormap(s2, 'hot');
%         colorbar
%         truesize(gca,[200 200]);
%         writeVideo(v,getframe(gca));
%     end
% end

activities = zeros(size(grid));
activities(grid > 0) = Max_Act;

for j = 173:1825
    j
    tic
    tailbudI = find(grid ~= borderValue);
    randomCells = randi([1 numel(tailbudI)],[1 numel(tailbudI)]);
    randomCells = tailbudI(randomCells);
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
        if ( grid(x(i), y(i) ) ~= rand_neighbor && grid(x(i),y(i)) ~= -1 && rand_neighbor ~= -1 )
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
            newH_J = H_J(newgrid, J_Cell_Cell, J_Cell_Medium, x(i), y(i));
            oldH_J = H_J(grid, J_Cell_Cell, J_Cell_Medium, x(i), y(i));
            newHp = sum(SA_constraint*(newp(allIs)-L0).^2);
            newHv = sum(volume_constraint*(newv(allIs)-A0).^2);
            newH = newHp+newHv;
            deltaH2 = newH - H - newH_J + oldH_J;
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
        g = grid;
        g(g < 0) = 0;
        imshow(g/max(g(:)));
        colormap(s1, 'colorcube');
        s2 = subplot(1,2,2);
        imshow(activities,[0 Max_Act]);
        colormap(s2, 'hot');
        colorbar
        truesize(gca,[200 200]);
        writeVideo(v,getframe(gca));
    end
    
    activities(activities ~= 0) = activities(activities ~= 0) - 1;
    if ( mod(j, 50) == 0 )
        grid = padarray(grid, 13, 0, 'post');
        tb = grid;
        tb(292:end,92) = borderValue;
        tb(292:end,111) = borderValue;
        tb(292:end,292) = borderValue;
        tb(292:end,311) = borderValue;
        tb(292,92:111) = borderValue;
        tb(292, 292:311) = borderValue;
        tb(292:end,92:111) = borderValue;
        tb(292:end, 292:311) = borderValue;
        grid = tb;
        currcell = max(grid(:)) + 1;
        for k = 1:14
            grid(end-12:end,111+12*(k-1):111+12*k - 1) = currcell;
            currcell = currcell + 1;
        end
        
        p = perim(grid,[]);
        p = p(p ~= -1);
        newp = p;
        volume = vol(grid);
        newv = volume;
        
        activities = padarray(activities, 13, 0, 'post');
        activities(end-12:end,111:278) = Max_Act;
    end
    toc
end
if (videoWrite)
    close(v);
end