videoWrite = 1;

rng(1);

if (videoWrite)
    v = VideoWriter('C:\Users\MH\Desktop\Notochord_Scoliosis.avi','Motion JPEG AVI');
    open(v);
end


grid = zeros(75,400);
cell_size = 20;

center_line = 100;

vacuolated = 1:30;
vacuolated(randi(30,1,10)) = [];

for i = 1:15
    grid(18:38,10+(i-1)*cell_size:10+i*cell_size) = i;
    grid(39:59,25+(i-1)*cell_size:25+i*cell_size) = i+15;
end

J_Cell_Medium = 20;
J_Cell_Cell = 65;
volume_constraint = 50;
SA_constraint = 25;
Temperature = 20;


A0_s = 325;
L0_s = 2*pi*sqrt(A0_s/pi)*8/pi;

p = perim_mex(grid,[]);
p(1) = [];
newp = p;
volume = vol(grid);
newv = volume;
gca = figure(1);

celltypes = zeros(size(grid,1),size(grid,2));
%Type 1 Cell Membrane
%Type 2 Cell Cytosol
%Type 3 Cell Vesicle
%Type 4 Cell Vacuole

celltypes(grid ~= 0) = 2;

membranes = zeros(size(grid,1),size(grid,2));
cells = unique(grid(:));

num_entities = max(grid(:));
A0 = repmat(A0_s,max(grid(:)),1);
L0 = repmat(L0_s,max(grid(:)),1);

for i = 1:length(cells)
    celltypes(bwperim(grid == i, 8)) = 1;
end

p_vacuole_to_vesicle = 0.05;
p_membrane_to_vesicle = 0.075;
p_vesicles_to_vacuole = 0.7;

for j = 1:50
    j
    randomCells = randi([1 numel(grid)],[1 numel(grid)]);
    [x,y] = ind2sub(size(grid),randomCells); 

    for i = 1:length(randomCells)
        [nx, ny] = meshgrid(x(i)-1:x(i)+1,y(i)-1:y(i)+1);
        n = [nx(:) ny(:)];
        n(n(:,1) == 0 | n(:,1) == size(grid,1)+1 | n(:,2) == 0 | n(:,2) == size(grid,2)+1,:) = [];
        inds = sub2ind(size(grid),n(:,1),n(:,2));
        all_neighbors = grid(inds);
        ri = randi([1 length(all_neighbors)]);
        rand_neighbor = all_neighbors(ri);
        ii = inds(ri);
        if ( grid(x(i), y(i) ) ~= 0 || rand_neighbor ~= 0)
            allIs = [grid(x(i), y(i) ) rand_neighbor];
            allIs(allIs == 0) = [];
            Hp = sum((p(allIs)-L0(allIs)).^2*SA_constraint);
            Hv = sum((volume(allIs)-A0(allIs)).^2*volume_constraint);
            H = Hp+Hv;
            newgrid = grid;
            newgrid(x(i), y(i) ) = rand_neighbor;
            prev = grid(x(i), y(i) );
            
            np = perim_mex(newgrid,[prev;rand_neighbor]);
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
            newHp = sum(SA_constraint*(newp(allIs)-L0(allIs)).^2);
            newHv = sum(volume_constraint*(newv(allIs)-A0(allIs)).^2);
            newH = newHp+newHv;
            deltaH2 = newH - H - newH_J_mex + oldH_J_mex;
            pr2 = exp(-deltaH2/Temperature);
            if ( pr2 > rand())             
                if (celltypes(inds(ri)) == 1 && prev == rand_neighbor && celltypes(x(i),y(i)) == 2 && ismember(prev, vacuolated))
                    if (rand() < p_membrane_to_vesicle)
                        celltypes(x(i),y(i)) = 3;
                    end
                else
                    if (~(celltypes(inds(ri)) == 2 && celltypes(x(i),y(i)) == 1) && celltypes(inds(ri)) ~= 1 && celltypes(x(i),y(i)) ~= 1 && grid(x(i),y(i)) == grid(inds(ri)))
                        celltypes(x(i),y(i)) = celltypes(inds(ri));
                    end

                    if (prev ~= 0)
                        previnds = find(grid == prev);
                        prevboundaries = find(bwperim(grid == prev));
                        for m = 1:length(previnds)
                            if (ismember(previnds(m), prevboundaries))
                                celltypes(previnds(m)) = -1;
                            end
                        end

                        celltypes(celltypes == 1 & grid == prev) = 2;
                        celltypes(grid == 0) = 0;
                        celltypes(celltypes == -1) = 1;
                        
                    end

                    if (rand_neighbor ~= 0)
                        neighinds = find(grid == rand_neighbor);
                        neighboundaries = find(bwperim(grid == rand_neighbor));
                        for m = 1:length(neighinds)
                            if (ismember(neighinds(m), neighboundaries))
                                celltypes(neighinds(m)) = -1;
                            end
                        end

                        celltypes(celltypes == 1 & grid == rand_neighbor) = 2;
                        celltypes(grid == 0) = 0;
                        celltypes(celltypes == -1) = 1;
                        
                    end

                    grid = newgrid;
                    p = newp;
                    volume = newv;
                end
            end
        end
    end
    
    for i = 1:length(unique(grid(:)))
        vesicles_pos = find(grid == i & celltypes == 3);
        for k = 1:length(vesicles_pos)
            [vesicle_pos_x,vesicle_pos_y] = ind2sub(size(grid), vesicles_pos(k));
            vesicle_pos_xy = [vesicle_pos_x, vesicle_pos_y];
            vesicle_pos_candidates = [vesicle_pos_xy(1) - 1 vesicle_pos_xy(2); vesicle_pos_xy(1) + 1 vesicle_pos_xy(2);...
                vesicle_pos_xy(1) vesicle_pos_xy(2) - 1; vesicle_pos_xy(1) vesicle_pos_xy(2) + 1];
            if (celltypes(vesicle_pos_candidates(1,1),vesicle_pos_candidates(1,2)) == 3 || celltypes(vesicle_pos_candidates(1,1),vesicle_pos_candidates(1,2)) == 4)
                if (rand() < p_vesicles_to_vacuole )
                    celltypes(vesicle_pos_candidates(1,1),vesicle_pos_candidates(1,2)) = 4;
                    celltypes(vesicle_pos_xy(1),vesicle_pos_xy(2)) = 4;
                end
            end
            if (celltypes(vesicle_pos_candidates(2,1),vesicle_pos_candidates(2,2)) == 3 || celltypes(vesicle_pos_candidates(2,1),vesicle_pos_candidates(2,2)) == 4)
                if (rand() < p_vesicles_to_vacuole )
                    celltypes(vesicle_pos_candidates(2,1),vesicle_pos_candidates(2,2)) = 4;
                    celltypes(vesicle_pos_xy(1),vesicle_pos_xy(2)) = 4;
                end
            end
            if (celltypes(vesicle_pos_candidates(3,1),vesicle_pos_candidates(3,2)) == 3 || celltypes(vesicle_pos_candidates(3,1),vesicle_pos_candidates(3,2)) == 4)
                if (rand() < p_vesicles_to_vacuole )
                    celltypes(vesicle_pos_candidates(3,1),vesicle_pos_candidates(3,2)) = 4;
                    celltypes(vesicle_pos_xy(1),vesicle_pos_xy(2)) = 4;
                end
            end
            if (celltypes(vesicle_pos_candidates(4,1),vesicle_pos_candidates(4,2)) == 3 || celltypes(vesicle_pos_candidates(4,1),vesicle_pos_candidates(4,2)) == 4)
                if (rand() < p_vesicles_to_vacuole )
                    celltypes(vesicle_pos_candidates(4,1),vesicle_pos_candidates(4,2)) = 4;
                    celltypes(vesicle_pos_xy(1),vesicle_pos_xy(2)) = 4;
                end
            end
        end
    end
    
    for i = 1:length(unique(grid(:)))
        vacuole_pos = find(bwperim(grid == i & celltypes == 4));
        for k = 1:length(vacuole_pos)
            if (rand() < p_vacuole_to_vesicle)
                [vacuole_pos_x,vacuole_pos_y] = ind2sub(size(grid), vacuole_pos(k));
                vacuole_pos_xy = [vacuole_pos_x,vacuole_pos_y];
                vacuole_pos_candidates = [vacuole_pos_xy(1) - 1 vacuole_pos_xy(2); vacuole_pos_xy(1) + 1 vacuole_pos_xy(2);...
                    vacuole_pos_xy(1) vacuole_pos_xy(2) - 1; vacuole_pos_xy(1) vacuole_pos_xy(2) + 1];
                
                vacuole_pos_candidates_inds = sub2ind(size(grid),vacuole_pos_candidates(:,1),vacuole_pos_candidates(:,2));
                vacuole_pos_candidates(celltypes(vacuole_pos_candidates_inds) ~= 2,:) = [];
                
                if (~isempty(vacuole_pos_candidates))
                    index = randi(size(vacuole_pos_candidates,1),1,1);
                    newpos = vacuole_pos_candidates(index,:);
                    oldcell = celltypes(newpos(1),newpos(2));
                    celltypes(newpos(1),newpos(2)) = 3;
                    celltypes(vacuole_pos_xy(1),vacuole_pos_xy(2)) = 2;
                end
            end
        end
    end
    
    for i = 1:length(unique(grid(:)))
        vesicles_pos = find(grid == i & celltypes == 3);
        for k = 1:length(vesicles_pos)
            [vesicle_pos_x,vesicle_pos_y] = ind2sub(size(grid), vesicles_pos(k));
            vesicle_pos_xy = [vesicle_pos_x, vesicle_pos_y];
            vesicle_pos_candidates = [vesicle_pos_xy(1) - 1 vesicle_pos_xy(2); vesicle_pos_xy(1) + 1 vesicle_pos_xy(2);...
                vesicle_pos_xy(1) vesicle_pos_xy(2) - 1; vesicle_pos_xy(1) vesicle_pos_xy(2) + 1];
            index = randi(4,1,1);
            newpos = vesicle_pos_candidates(index,:);
            oldcell = celltypes(newpos(1),newpos(2));
            if (oldcell ~= 1)
                celltypes(newpos(1),newpos(2)) = 3;
                celltypes(vesicle_pos_xy(1),vesicle_pos_xy(2)) = oldcell;
            else
                celltypes(vesicle_pos_xy(1),vesicle_pos_xy(2)) = 2;
            end
        end
    end
    
    for i = 1:length(unique(grid))
        A0(i) = A0_s + sum(sum(grid == i & (celltypes == 3 | celltypes == 4)));
        L0(i) = 2*pi*sqrt(A0(i)/pi);
    end
    
    gca = figure(1);
    imshow(grid/max(grid(:)));
    colormap('colorcube');
    if (videoWrite)
        writeVideo(v,getframe(gca));
    end
end
% 
% for j = 1:300
%     j
%     randomCells = randi([1 numel(grid)],[1 numel(grid)]);
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
%         if ( grid(x(i), y(i) ) ~= rand_neighbor )
%             allIs = [grid(x(i), y(i) ) rand_neighbor];
%             allIs(allIs == 0) = [];
%             Hp = sum((p(allIs)-L0).^2*SA_constraint);
%             Hv = sum((volume(allIs)-A0).^2*volume_constraint);
%             H = Hp+Hv;
%             newgrid = grid;
%             newgrid(x(i), y(i) ) = rand_neighbor;
%             prev = grid(x(i), y(i) );
% 
%             newp = p;
%             np = perim_mex(newgrid,[prev rand_neighbor]);
%             if (np(1) ~= -1)
%                 newp(prev) = np(1);
%             end
%             if (np(2) ~= -1)
%                 newp(rand_neighbor) = np(2);
%             end
% 
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
%             newH_J_mex = H_J_mex(newgrid, J_Cell_Cell, J_Cell_Medium, x(i), y(i));
%             oldH_J_mex = H_J_mex(grid, J_Cell_Cell, J_Cell_Medium, x(i), y(i));
%             newHp = sum(SA_constraint*(newp(allIs)-L0).^2);
%             newHv = sum(volume_constraint*(newv(allIs)-A0).^2);
%             newH = newHp+newHv;
%             deltaH2 = newH - H - newH_J_mex + oldH_J_mex;
%             pr2 = exp(-deltaH2/Temperature);
%             if ( pr2 > rand() )
%                 grid = newgrid;
%                 p = newp;
%                 volume = newv;
%                 H = newH;
%             end
%         end
%     end
%     
%     if (videoWrite)
%         gca = figure(1);
%         imshow(grid/max(grid(:)));
%         colormap('colorcube');
%         writeVideo(v,getframe(gca));
%     end
% end
% if (videoWrite)
%     close(v);
% end
% imshow(grid);