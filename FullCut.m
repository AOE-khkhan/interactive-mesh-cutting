clear all; close all;
finCharName = 'CaptainAmerica';
charName = finCharName;
th = 8;
% variables = load(['TensorField_', charName]);
% load(['TensorField_', charName]);
load(['ProcessedAll_', charName]);
charName = finCharName;
load(['MeshThinning_', charName, '_th', sprintf('%02d', th)]);
load(['MeshConnections_', charName]);

charName = [charName, '_20161024'];
prompt = 'Enter your name: ';
cutterName = [input(prompt, 's')];
fName = [charName, '_by_', cutterName];
colLabels = jet(20); close all; 
vertices_orig = vertices;
views = [0 90; 180 270];
thCond = 1;

%% View different pathes at specific threshold
for subp = 1:2
    figure(subp); hold on;
    trisurf(newfaces, vertices(:,1), vertices(:,2), vertices(:,3), 'FaceColor', [0.95,0.95,1.0], 'EdgeColor', [0.1,0.1,0.1]);
    view(views(subp, :)); %hold on;
    pbaspect([1 2 1]);
    set(gcf, 'units', 'normalized', 'outerposition', [(subp-1)*0.5 0 0.5 1]);
end

for clus = 1:length(patchEdgesPerClus)
    patchedges = patchEdgesPerClus{clus}                             %% Saves the patch edges per cluster
    for subp = 1:2
        figure(subp);hold on;
        for i = 1:size(patchedges,1),
            curredge = patchedges(i,:);
            x = line(vertices([curredge(1) curredge(2)],1),vertices([curredge(1) curredge(2)],2), vertices([curredge(1) curredge(2)],3));
            set(x, 'LineWidth', 5);
            set(x, 'Color', 'b');
        end
    end
end

%% Re-adjusts the values of D2 (curvature) for shortest path algorithm
D2_orig = D2;
D2_weighted = abs(-(min(D2)-D2));
% D2 = -(min(D2)-D2);

for i = 1:size(patchedges,1),
   curredge = patchedges(i,:);
   x = line(vertices([curredge(1) curredge(2)],1),vertices([curredge(1) curredge(2)],2), vertices([curredge(1) curredge(2)],3))
   set(x, 'Color', colLabels(th,:));
end

% for subp = 3:4
%     figure(subp);
%     trisurf(newfaces, vertices(:,1), vertices(:,2), vertices(:,3), 'FaceColor', [0.95,0.95,1.0], 'EdgeColor', [0.1,0.1,0.1]);
%     view(views(subp-2, :)); hold on;
%     for clus = 1:size(patchEdgesPerClus,2)
% %         patchedgescopy = patchedges;
%         patchedges = patchEdgesPerClus{clus};
%         for i = 1:size(patchedges,1),
%            curredge = patchedges(i,:);
%            x = line(vertices([curredge(1) curredge(2)],1),vertices([curredge(1) curredge(2)],2), vertices([curredge(1) curredge(2)],3))
%            set(x, 'Color', 'b', 'Linewidth', 2);
%         end
%     end
%     figure(subp); pbaspect([1 2 1]);
%     set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]);
% end

%% Mesh cut selection
numCuts = 0;
finalPaths = []; figNum = 1; indexRev = [];
finalIndex = [];
prompt1 = 'Full cut or slice? [1/0]: ';
prompt = 'Number of cuts: ';
meshCut = input(prompt, 's');
meshCut = str2num(meshCut);
typeOfCut = [];
close all;
for cutNum = 1:meshCut
    typeOfCut{cutNum} = input(['Cut ', num2str(cutNum), ' - ', prompt1], 's');
end

for cutNum = 1:meshCut
    numCuts = numCuts + 1;
%     pause
%     typeOfCut{cutNum} = input(prompt1, 's');

    for slices = 1:2
%         figNum = figNum + 1;
        figure(slices); hold on;
        trisurf(newfaces, vertices(:,1), vertices(:,2), vertices(:,3), 'FaceColor', ...
            [0.95,0.95,1.0], 'EdgeColor', [0.1,0.1,0.1], 'FaceAlpha', 0.7);
        view(views(slices, :)); %hold on;
        indexSel = finalIndex(:); indexSel = indexSel(indexSel~=0);
        scatter3(vertices(indexSel,1),vertices(indexSel,2),vertices(indexSel,3), 150, 'o', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');
        
        for clus = 1:size(patchEdgesPerClus,2)
    %         patchedgescopy = patchedges;
            patchedges = patchEdgesPerClus{clus};
            for i = 1:size(patchedges,1),
               curredge = patchedges(i,:);
               x = line(vertices([curredge(1) curredge(2)],1),vertices([curredge(1) curredge(2)],2), vertices([curredge(1) curredge(2)],3));
               set(x, 'Color', 'b', 'Linewidth', 2);
            end
        end
        figure(slices); pbaspect([1 2 1]);
        set(gcf, 'units', 'normalized', 'outerposition', [(slices-1)*0.5 0 0.5 1]);
         
        disp('Zoom in or rotate on region of interest, then press Return.');
        pause
        dcm_obj = datacursormode(figure(slices)); disp('Select two (2) vertices...');
        set(dcm_obj,'DisplayStyle','datatip', 'SnapToDataVertex','on','Enable','on')
        disp('Click line to display a data tip, then press Return.');
        pause
        
%         tic
        
        c_info = getCursorInfo(dcm_obj);
        index = [];
        for i = 1:size(c_info,2),
            curpoint = c_info(i).Position;
            vertices1 = vertices; %(indix2,:);
            currindex1 = find(vertices1(:,1) == curpoint(1));
            currindex2 = find(vertices1(currindex1,2) == curpoint(2));
            currindex3 = find(vertices1(currindex1(currindex2),3) == curpoint(3));
            index(i) = currindex1(currindex2(currindex3));
        end
        indexRev = index(end:-1:1);
        finalIndex(numCuts,slices,1:2) = indexRev;
        [dists,paths] = dijkstra1(nodesV,segsV,D2_weighted, indexRev(1), indexRev(2));
        finalPaths{numCuts, slices} = paths;
        close all;
%         toc
    end
    [dists,paths] = dijkstra1(nodesV,segsV,D2_weighted,finalIndex(numCuts, 1,2), finalIndex(numCuts, 2,2));
    finalPaths{numCuts, 3} = paths;
    
    if typeOfCut{cutNum} == '1'
        [dists,paths] = dijkstra1(nodesV,segsV,D2_weighted,finalIndex(numCuts, 2,1), finalIndex(numCuts, 1,1));
        finalPaths{numCuts, 4} = paths;
    end
    close all;
end
allPaths = finalPaths; 

figure(5); hold on;
trisurf(newfaces, vertices(:,1), vertices(:,2), vertices(:,3), 'FaceColor', [0.95,0.95,1.0], 'EdgeColor', [0.1,0.1,0.1]);
for i=1:size(allPaths,1)
    for j=1:3
        index = allPaths{i,j};
        scatter3(vertices(index,1),vertices(index,2),vertices(index,3), 100, 'o', 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'g');
    end
end
pbaspect([1 2 1]);

for i = 1:2
    figure(5); pbaspect([1 2 1]); view(views(i,:));
    fname = [charName, '_', cutterName, '_th', sprintf('%02d', th), '_Cut_View', num2str(i)];
    set(gcf, 'units', 'normalized', 'outerposition', [(i-1)*0.5 0 0.5 1]);
    printPNG(figure(5), [fname, '.png']);
end
fname = [charName, '_', cutterName, '_th', sprintf('%02d', th), '_Cut_Views'];
printPNG(figure(5), [fname, '.fig']);

tic

%% Deletes path vertices for cutting
% finalPaths = [];
% numPath = 1;
% for i=1:size(allPaths, 2)
%     finalPaths{i} = allPaths{numPath,i};
% end
% finalPath1 = finalPaths{1};
% finalPath2 = finalPaths{2};
% finalPath3 = finalPaths{3};
% finPaths = [finalPath1 finalPath3(2:end-1) finalPath2(end:-1:1) ];
% 

newfaces1 = newfaces;
newedges1 = newedges;
newfaces2 = newfaces;
for cutnum = 1:numCuts,
    for i = 1:size(allPaths,2),
        paths = allPaths{cutnum, i};
        for j =  1:size(paths,2),
            currvert = paths(j);
            currvertindex = find(newfaces2 == currvert);
            newfaces2(mod(currvertindex, size(newfaces2,1)),:) = [];
        end
    end
end
% 
% edgecounter = 1;
% for i = 1:size(newfaces1,1),
%     edges1(edgecounter,:) = sort([newfaces1(i,1) newfaces1(i,2)]);
%     oppositevertex1(edgecounter,1) = newfaces1(i,3);
%     faceref1(edgecounter,1) = i;
%     edgecounter = edgecounter + 1;
%     edges1(edgecounter,:) = sort([newfaces1(i,2) newfaces1(i,3)]);
%     oppositevertex1(edgecounter,1) = newfaces1(i,1);
%     faceref1(edgecounter,1) = i;
%     edgecounter = edgecounter + 1;
%     edges1(edgecounter,:) = sort([newfaces1(i,1) newfaces1(i,3)]);
%     oppositevertex1(edgecounter,1) = newfaces1(i,2);
%     faceref1(edgecounter,1) = i;
%     edgecounter = edgecounter + 1;
%     
%     
%     face = newfaces1(i,:);
%     facenormals1(i,:) = cross(vertices(face(2),:)-vertices(face(1),:), vertices(face(3),:)-vertices(face(1),:));
%     facearea1(i,1) = (1/2)*norm(facenormals1(i,:));
% end
% 
% edgecounter = 1;
% newedges1 = [0 0];
% newfaceref1 = [0 0];
% newoppositevertex1 = [0 0];
% for i = 1:size(edges1,1),
%     edge = edges1(i,:);
%     indexfound = [];
%     if i ~= 1,
%         indexrows = find(newedges1(:,1) == edge(1));
%         indexfound = newedges1(newedges1(indexrows,2) == edge(2));
%         index = find(newedges1(indexrows,2) == edge(2));
%         index = indexrows(index);
%     end
%     if length(indexfound) >= 1,
%         newfaceref1(index,2) = faceref1(i);
%         newoppositevertex1(index,2) = oppositevertex1(i);
%         disp('FOUND!');
%     else
%         disp('NOT FOUND!');
%         newedges1(edgecounter,:) = edge;
%         newfaceref1(edgecounter,1) = faceref1(i);
%         newoppositevertex1(edgecounter,1) = oppositevertex1(i);
%         edgecounter = edgecounter + 1;
%     end
% end
% 
% edgeangle1 = [];
% 
% for i = 1:size(newedges1,1),
%     twofaces = newfaceref1(i,:);
%     if (twofaces(2) ~= 0),
%         normal1 = facenormals1(twofaces(1),:);
%         normal2 = facenormals1(twofaces(2),:);
%         edgeangle1(i) = atand(norm(cross(normal1,normal2))/(dot(normal1,normal2)));
%     end
%     newedgeslength1(i,1) = norm(vertices(newedges1(i,1), :) - vertices(newedges1(i,2), :)); 
% end

%% Mesh Cutting 
newfacescopy = newfaces;
newedgescopy = newedges;

vertices2 = vertices;
connectingFaces = [];
finalfaces = newfaces2;
for cutnum = 1:numCuts
    finalPaths = [];
    numPath = 1;
    for i=1:size(allPaths,2)
        finalPaths{i} = allPaths{cutnum,i};
    end
    finalPath1 = finalPaths{1};
    finalPath2 = finalPaths{2};
    finalPath3 = finalPaths{3};
    if typeOfCut{cutnum} == '1'
        finalPath4 = finalPaths{4};
        finPaths = [finalPath1 finalPath3(2:end-1) finalPath2(end:-1:1) finalPath4(2:end)];
    else
        finPaths = [finalPath1 finalPath3(2:end-1) finalPath2(end:-1:1) ];
    end
    
    %% Delete cut path from faces
    newfaces = newfacescopy;
    newedges = newedgescopy;
    newfaces1 = newfacescopy;
    %newedges1 = newedgescopy;
    %newoppositevertex1 = newoppositevertex;
    for i = 1:size(allPaths,2)
        paths = allPaths{cutnum, i};
        for j =  1:size(paths,2),
            currvert = paths(j);
            currvertindex = find(newfaces1 == currvert);
            newfaces1(mod(currvertindex, size(newfaces1,1)),:) = [];
        end
    end
    
    
    edgecounter = 1;
    for i = 1:size(newfaces1,1),
        edges1(edgecounter,:) = sort([newfaces1(i,1) newfaces1(i,2)]);
        oppositevertex1(edgecounter,1) = newfaces1(i,3);
        faceref1(edgecounter,1) = i;
        edgecounter = edgecounter + 1;
        edges1(edgecounter,:) = sort([newfaces1(i,2) newfaces1(i,3)]);
        oppositevertex1(edgecounter,1) = newfaces1(i,1);
        faceref1(edgecounter,1) = i;
        edgecounter = edgecounter + 1;
        edges1(edgecounter,:) = sort([newfaces1(i,1) newfaces1(i,3)]);
        oppositevertex1(edgecounter,1) = newfaces1(i,2);
        faceref1(edgecounter,1) = i;
        edgecounter = edgecounter + 1;
        
        
        face = newfaces1(i,:);
        facenormals1(i,:) = cross(vertices(face(2),:)-vertices(face(1),:), vertices(face(3),:)-vertices(face(1),:));
        facearea1(i,1) = (1/2)*norm(facenormals1(i,:));
    end
    
    edgecounter = 1;
    newedges1 = [0 0];
    newfaceref1 = [0 0];
    newoppositevertex1 = [0 0];
    for i = 1:size(edges1,1),
        edge = edges1(i,:);
        indexfound = [];
        if i ~= 1,
            indexrows = find(newedges1(:,1) == edge(1));
            indexfound = newedges1(newedges1(indexrows,2) == edge(2));
            index = find(newedges1(indexrows,2) == edge(2));
            index = indexrows(index);
        end
        if length(indexfound) >= 1,
            newfaceref1(index,2) = faceref1(i);
            newoppositevertex1(index,2) = oppositevertex1(i);
            disp('FOUND!');
        else
            disp('NOT FOUND!');
            newedges1(edgecounter,:) = edge;
            newfaceref1(edgecounter,1) = faceref1(i);
            newoppositevertex1(edgecounter,1) = oppositevertex1(i);
            edgecounter = edgecounter + 1;
        end
    end
    
    edgeangle1 = [];
    
    for i = 1:size(newedges1,1),
        twofaces = newfaceref1(i,:);
        if (twofaces(2) ~= 0),
            normal1 = facenormals1(twofaces(1),:);
            normal2 = facenormals1(twofaces(2),:);
            edgeangle1(i) = atand(norm(cross(normal1,normal2))/(dot(normal1,normal2)));
        end
        newedgeslength1(i,1) = norm(vertices(newedges1(i,1), :) - vertices(newedges1(i,2), :));
    end

    paths = finPaths;
    firstedge = sort(paths(1:2));
    edgeindex1 = find(newedges(:,1) == firstedge(1));
    edgeindex2 = find(newedges(edgeindex1,2) == firstedge(2));
    opvertsfirst = newoppositevertexcopy(edgeindex1(edgeindex2),:);
    
    opvert = opvertsfirst(1);
    
    firstedge = sort(paths(end-1:end));
    edgeindex1 = find(newedgescopy(:,1) == firstedge(1));
    edgeindex2 = find(newedgescopy(edgeindex1,2) == firstedge(2));
    opvertsend = newoppositevertexcopy(edgeindex1(edgeindex2),:);
    
    index = find(newfaceref1(:,2) == 0);
    holeedges = newedges1(index,:);
    
    startindex = find(holeedges(:,1) == opvertsfirst(1));
    edgecopy = holeedges(startindex,:);
    holeedges(startindex,:) = [];
    holeedges = [edgecopy; holeedges];
    
    if typeOfCut{cutnum} == '0',
        
        counter = 1;
        sortholeedges = [];
        while (size(holeedges,1) ~= 0),
            if (counter == 1),
                edge = holeedges(1,:);
                holeedges(1,:) = [];
            else
                edge = sortholeedges(counter,:);
            end
            sortholeedges(counter,:) = edge;
            counter = counter + 1;
            endedge = edge(2);
            endedgefound = find(holeedges == endedge);
            if (size(endedgefound,1) ~= 1),
                endedgefound = endedgefound(1);
            end
            if endedgefound <= size(holeedges,1),
                sortholeedges(counter,:) = holeedges(endedgefound,:);
                holeedges(endedgefound,:) = [];
            else
                sortholeedges(counter,:) = holeedges(endedgefound-size(holeedges,1),end:-1:1);
                holeedges(endedgefound-size(holeedges,1),:) = [];
            end
        end
        endindex = find(sortholeedges(:,2) == opvertsend(1)); % opvertsend(2) before
        firsthole = sortholeedges(1:endindex,:);
        secondhole = sortholeedges(endindex+1:end,:);
        connectingFaces = [connectingFaces; [firsthole(1,1), secondhole(end,1), secondhole(1,1)]];
        connectingFaces = [connectingFaces; [firsthole(1,1), firsthole(end,1), secondhole(end,1)]];
        pairs1 = firsthole(:,1);
        pairs2 = secondhole(:,1);
    else
        counter = 1;
        sortholeedges = [];
        while (size(holeedges,1) ~= 0),
            if (counter == 1),
                edge = holeedges(1,:);
                holeedges(1,:) = [];
            else
                edge = sortholeedges(counter,:);
            end
            sortholeedges(counter,:) = edge;
            counter = counter + 1;
            endedge = edge(2);
            endedgefound = find(holeedges == endedge);
            if (size(endedgefound,1) == 0),
                firsthole = sortholeedges;
               break; 
            end
            if (size(endedgefound,1) ~= 1),
                endedgefound = endedgefound(1);
            end
            if endedgefound <= size(holeedges,1),
                sortholeedges(counter,:) = holeedges(endedgefound,:);
                holeedges(endedgefound,:) = [];
            else
                sortholeedges(counter,:) = holeedges(endedgefound-size(holeedges,1),end:-1:1);
                holeedges(endedgefound-size(holeedges,1),:) = [];
            end
        end
        counter = 1;
        sortholeedges = [];
        while (size(holeedges,1) ~= 0),
            if (counter == 1),
                edge = holeedges(1,:);
                holeedges(1,:) = [];
            else
                edge = sortholeedges(counter,:);
            end
            sortholeedges(counter,:) = edge;
            counter = counter + 1;
            endedge = edge(2);
            endedgefound = find(holeedges == endedge);
            if (size(endedgefound,1) == 0),
               secondhole = sortholeedges;
               break; 
            end
            if (size(endedgefound,1) ~= 1),
                endedgefound = endedgefound(1);
            end
            if endedgefound <= size(holeedges,1),
                sortholeedges(counter,:) = holeedges(endedgefound,:);
                holeedges(endedgefound,:) = [];
            else
                sortholeedges(counter,:) = holeedges(endedgefound-size(holeedges,1),end:-1:1);
                holeedges(endedgefound-size(holeedges,1),:) = [];
            end
        end
        secondhole = sortholeedges;
        pairs1 = firsthole(:,1);
        pairs2 = secondhole(:,1);
    end
    %% Fill hole
    for l = 1:2,
        if l == 1,
        sortholeedges = pairs1;
        else
        sortholeedges = pairs2;
        end

        vertices1 = vertices([sortholeedges(:,1); sortholeedges(1)], :);
        W = [];
        mu = [];
        for i = 1:size(vertices1,1)-1,
            W(i,i+1) = 0;
            mu(i,i+1) = 0;
        end

        for i = 1:size(vertices1,1)-2,
            v1 = vertices1(i,:) - vertices1(i+1,:);
            v2 = vertices1(i+2,:) - vertices1(i+1,:);
            vcross = cross(v1,v2);
            W(i, i + 2) = 1/2*norm(vcross);
            neighbors1 = find(newedges1(:,1) == sortholeedges(i+1,1));
            neighbors2 = find(newedges1(:,2) == sortholeedges(i+1,1));
            mu(i, i+2) = max(edgeangle1([neighbors1; neighbors2]));
        end

        j = 2;

        O = zeros(size(vertices1,1),size(vertices1,1));

        while j < size(vertices1,1),
            j = j+1;
            for i = 1: size(vertices1,1)-j,
                k = i+j;
                counter = 1;
                M = [];
                setW = [];
                setmu =[];
                for m = i+1:k-1,

                    dihed = 0;
                    v1 = vertices1(i,:) - vertices1(m,:);
                    v2 = vertices1(k,:) - vertices1(m,:);
                    vcross = cross(v1,v2);

                    if (O(i,m)~=0 && O(m,k)~=0),
                        v3 = vertices1(i,:) - vertices1(O(i,m),:);
                        v4 = vertices1(m,:) - vertices1(O(i,m),:);
                        cross1 = cross(v3,v4);


                        v5 = vertices1(m,:) - vertices1(O(m,k),:);
                        v6 = vertices1(k,:) - vertices1(O(m,k),:);
                        cross2 = cross(v5,v6);

                        dihed = atand(norm(cross(cross1,cross2))/(dot(cross1,cross2)));
                    end

                    setmu(counter) = max([mu(i,m) mu(m,k) dihed]);
                    setW(counter) = W(i,m) + W(m,k) + 1/2*norm(vcross);
                    M(counter) = m;
                    counter = counter + 1;
                end

                mu(i,k) = min(setmu);
                indix = find(setmu == min(setmu));
                if (size(indix,2) ~= 1),
                    setW1 = setW(indix);
                    indix1 = find(setW1 == min(setW1));
                    indix1 = max(indix1);
                    W(i,k) = setW(indix(indix1));
                    O(i,k) = M(indix(indix1));
                else
                    W(i,k) = setW(max(indix));
                    O(i,k) = M(indix(find(indix == max(indix))));
                end
            end
        end

        S = [];
        S = Trace(1, size(vertices1,1), O, S);
        S1 = [];
        counter = 1;
        for i = 1:3:size(S,2)-2,
            S1(counter,:) = [S(i) S(i+1) S(i+2)];
            counter = counter + 1;
        end

        S2 = S1;
        %S2(find(S2 == size(sortholeedges,1))) = 1;
        S2(find(S2 == size(sortholeedges,1)+1)) = 1;

        Sr = S2;
        Sr = unique(Sr, 'rows');
        counter = 1;
        rowstodelete = [];
        for i = 1:size(Sr,1),
            row = Sr(i,:);
            if (length(unique(row)) ~= 3),
                rowstodelete(counter) = i;
                counter = counter + 1;
            end
        end

        Sr(rowstodelete,:) = [];
        
        Sr = sortholeedges(Sr);
        
        %% Refinement
        Srcopy = Sr;
        S21 = Sr;
        interioredgescounter = 1;
        for i = 1:size(S21,1),
            interioredges(interioredgescounter,:) = sort([S21(i,1) S21(i,2)]);
            interiorfaceref(interioredgescounter,1) = i;
            interioroppositevertex(interioredgescounter,1) = S21(i,3);
            interioredgescounter = interioredgescounter + 1;
            interioredges(interioredgescounter,:) = sort([S21(i,2) S21(i,3)]);
            interiorfaceref(interioredgescounter,1) = i;
            interioroppositevertex(interioredgescounter,1) = S21(i,1);
            interioredgescounter = interioredgescounter + 1;
            interioredges(interioredgescounter,:) = sort([S21(i,1) S21(i,3)]);
            interiorfaceref(interioredgescounter,1) = i;
            interioroppositevertex(interioredgescounter,1) = S21(i,2);
            interioredgescounter = interioredgescounter + 1;
        end
        
        
        
        for i = 1:size(vertices,1),
            neighbors1 = find(newedges(:,1) == i);
            neighbors2 = find(newedges(:,2) == i);
            aveedgelength(i,1) = mean(newedgeslength([neighbors1; neighbors2],1));
        end
        
        alphafactor = sqrt(2);
        newinterioredgecounter = 1;
        newinterioredges = [0 0];
        newinteriorfaceref = [0 0];
        newinterioroppositevertex = [0 0];
        for i = 1:size(interioredges,1),
            edge = interioredges(i,:);
            indexfound = [];
            if i ~= 1,
                indexrows = find(newinterioredges(:,1) == edge(1));
                indexfound = newinterioredges(newinterioredges(indexrows,2) == edge(2));
                index = find(newinterioredges(indexrows,2) == edge(2));
                index = indexrows(index);
            end
            if length(indexfound) >= 1,
                newinteriorfaceref(index,2) = interiorfaceref(i);
                newinterioroppositevertex(index,2) = interioroppositevertex(i);
                disp('FOUND!');
            else
                disp('NOT FOUND!');
                newinterioredges(newinterioredgecounter,:) = edge;
                newinteriorfaceref(newinterioredgecounter,1) = interiorfaceref(i);
                newinterioroppositevertex(newinterioredgecounter,1) = interioroppositevertex(i);
                newinterioredgecounter = newinterioredgecounter + 1;
            end
        end
        
        for i = 1:size(newinterioredges,1),
            newinterioredgeslength(i,1) = norm(vertices(newinterioredges(i,1),:) - vertices(newinterioredges(i,2),:));
        end
        
        
        counter = 1;
        rowstodelete = [];
        rcounter = 1;
        Sr = [];
%         newfacescopy = newfaces;
%         newedgescopy = newedges;
%         newedges = newedges1;
%         newfaces = newfaces1;
%         newedges1 = newedges;
%         newfaces1 = newfaces;
        newinterioredges1 = newinterioredges;
        %newoppositevertexcopy = newoppositevertex;
        newoppositevertex = newoppositevertex1;
        newoppositevertex1 = newoppositevertex;
        newinterioroppositevertex1 = newinterioroppositevertex;
        %sortholeedges1 = sortholeedges;
        %S21 = sortholeedges1(S2);
        
        %for i = 10:20,
        iterator = [1 2 3];
        %k = 1;
        while (size(S21,1) ~= 0),
            disp('i');
            i = 1;
            index = [S21(i,1) S21(i,2) S21(i,3)];
            disp(index);
            S21(i,:) = [];
            vc(i,:) = mean(vertices2(index,:));
            sigmavc = mean(aveedgelength(index));
            for k = 1:3,
                %k = 1
                if (alphafactor*(norm(vc(i,:) - vertices2(index(k),:))) > sigmavc) && (alphafactor*(norm(vc(i,:) - vertices2(index(k),:))) > aveedgelength(index(k))),
                    rowstodelete(counter) = i;
                    if (vertices2(end,:) ~= vc(i,:)),
                        disp('AHA');
                        vertices2(end+1,:) = vc(i,:);
                    end
                    % RELAX EDGES ij, jk, ik
                    edgerelaxed = 0;
                    for j = 1:3,
                        disp('j');
                        disp(j);
                        currentiterator = iterator(iterator ~= j);
                        currentedge = sort(index(currentiterator));
                        edgeindex1 = find(newinterioredges1(:,1) == currentedge(1));
                        edgeindex2 = find(newinterioredges1(edgeindex1,2) == currentedge(2));
                        if (size(edge,1) ~= 0),
                            outervertindex1 = find(newedges1(:,1) == edge(1));
                            outervertindex2 = find(newedges1(outervertindex1,2) == edge(2));
                            edge = newinterioredges1(edgeindex1(edgeindex2),:);
                        end
                        if (size(edge,1) == 0),
                            disp('EDGE NOT FOUND!!!!');
                            %disp([index2(j) index2(j+1)]);
                            break;
                        end
                        opverts = newinterioroppositevertex1(edgeindex1(edgeindex2),:);
                        if (opverts(2) == 0),
                            outervertindex1 = find(newedges1(:,1) == edge(1));
                            outervertindex2 = find(newedges1(outervertindex1,2) == edge(2));
                            opvert = newoppositevertex1(outervertindex1(outervertindex2),1);
                            %                     if (size(opvert,2) == 0),
                            %                         break;
                            %                     end
                            if (size(opvert,2) == 0),
                                break;
                            elseif (size(opvert,1) == 0),
                                break;
                            elseif (size(opvert,2) == 2),
                                break;
                            end
                        else
                            opvert = setdiff(index, edge);
                            opvert = opverts(opverts ~= opvert);
                            if (size(opvert,2) == 0),
                                break;
                            elseif (size(opvert,1) == 0),
                                break;
                            elseif (size(opvert,2) == 2),
                                break;
                            end
                        end
                        
                        edgelenij = norm(vertices2(edge(1),:) - vertices2(edge(2),:));
                        edgelenopp = norm(vc(i,:) - vertices2(opvert,:));
                        
                        midpointij = mean([vertices2(edge(1),:); vertices2(edge(2),:)],1);
                        
                        if (norm(vertices2(opvert,:) - midpointij) < edgelenij/2) && (norm(vertices2(size(vertices2,1),:) - midpointij) < edgelenij/2),
                            %if edgelenij > edgelenopp,
                            edgerelaxed = 1;
                            disp('EDGE RELAXED');
                            disp(edge)
                            disp('opvert');
                            disp(opverts);
                            disp(opvert);
                            if (opverts(2) == 0),
                                newedges1(outervertindex1(outervertindex2),:) = sort([opvert size(vertices2,1)]);
                                newoppositevertex1(outervertindex1(outervertindex2),:) = edge;
                            end
                            newinterioredges1(edgeindex1(edgeindex2),:) = sort([opvert size(vertices2,1)]);
                            newinterioroppositevertex1(edgeindex1(edgeindex2),:) = edge;
                            newinterioredgeslength(edgeindex1(edgeindex2),1) = norm(vertices2(edge(1),:) - vertices2(edge(2),:));
                            
                            if (opverts(2) == 0),
                                %newfaces1(newfaceref(outervertindex1(outervertindex2),1),:) = [opvert edge(1) size(vertices2,1)];
                                Sr(rcounter,:) = sort([opvert edge(1) size(vertices2,1)]);
                                rcounter = rcounter + 1;
                                Sr(rcounter,:) = sort([opvert edge(2) size(vertices2,1)]);
                                rcounter = rcounter + 1;
                            else
                                %Sr(rcounter-j,:) = [opvert index(j+1) size(vertices2,1)];
                                Sr(rcounter,:) = sort([opvert edge(1) size(vertices2,1)]);
                                rcounter = rcounter + 1;
                                Sr(rcounter,:) = sort([opvert edge(2) size(vertices2,1)]);
                                rcounter = rcounter + 1;
                            end
                            
                            tobedeletedface = sort([edge opvert]);
                            S21copy = sort(S21,2);
                            newfaces1copy = sort(newfaces1,2);
                            tobedeletedfaceindex1 = find(S21copy(:,1) == tobedeletedface(1));
                            tobedeletedfaceindex2 = find(S21copy(tobedeletedfaceindex1,2) == tobedeletedface(2));
                            tobedeletedfaceindex3 = find(S21copy(tobedeletedfaceindex1(tobedeletedfaceindex2),3) == tobedeletedface(3));
                            
                            
                            tobedeletednewfaceindex1 = find(newfaces1copy(:,1) == tobedeletedface(1));
                            tobedeletednewfaceindex2 = find(newfaces1copy(tobedeletednewfaceindex1,2) == tobedeletedface(2));
                            tobedeletednewfaceindex3 = find(newfaces1copy(tobedeletednewfaceindex1(tobedeletednewfaceindex2),3) == tobedeletedface(3));
                            
                            
                            
                            tobedeletedSrindex1 = find(Sr(:,1) == tobedeletedface(1));
                            tobedeletedSrindex2 = find(Sr(tobedeletedSrindex1,2) == tobedeletedface(2));
                            tobedeletedSrindex3 = find(Sr(tobedeletedSrindex1(tobedeletedSrindex2),3) == tobedeletedface(3));
                            
                            
                            
                            disp('tobedeleted');
                            disp(S21(tobedeletedfaceindex1(tobedeletedfaceindex2(tobedeletedfaceindex3)),:));
                            
                            disp('tobedeleted');
                            disp(Sr(tobedeletedSrindex1(tobedeletedSrindex2(tobedeletedSrindex3)),:));
                            
                            disp('tobedeleted');
                            disp(newfaces1(tobedeletednewfaceindex1(tobedeletednewfaceindex2(tobedeletednewfaceindex3)),:));
                            
                            S21(tobedeletedfaceindex1(tobedeletedfaceindex2(tobedeletedfaceindex3)),:) = [];
                            newfaces1(tobedeletednewfaceindex1(tobedeletednewfaceindex2(tobedeletednewfaceindex3)),:) = [];
                            Sr(tobedeletedSrindex1(tobedeletedSrindex2(tobedeletedSrindex3)),:) = [];
                            %S21 = [opvert edge(2) size(vertices2,1); S21];
                            
                            edge1 = sort([edge(1) opvert]);
                            edge2 = sort([edge(2) opvert]);
                            edge3 = sort([edge(1) size(vertices2,1)]);
                            edge4 = sort([edge(2) size(vertices2,1)]);
                            edgeset = [edge1; edge2; edge3; edge4];
                            for l = 1:size(edgeset,1),
                                curredge = edgeset(l,:);
                                refinededgeindex1 = find(newinterioredges1(:,1) == curredge(1));
                                refinededgeindex2 = find(newinterioredges1(refinededgeindex1,2) == curredge(2));
                                curroppvertex = newinterioroppositevertex1(refinededgeindex1(refinededgeindex2),:);
                                
                                if (l==1) || (l==3),
                                    currincludededgeindex = find(curroppvertex == edge(2));
                                    if (l==1),
                                        newinterioroppositevertex1(refinededgeindex1(refinededgeindex2), currincludededgeindex) = size(vertices2,1);
                                    else
                                        newinterioroppositevertex1(refinededgeindex1(refinededgeindex2), currincludededgeindex) = opvert;
                                    end
                                else
                                    currincludededgeindex = find(curroppvertex == edge(1));
                                    if (l==2),
                                        newinterioroppositevertex1(refinededgeindex1(refinededgeindex2), currincludededgeindex) = size(vertices2,1);
                                    else
                                        newinterioroppositevertex1(refinededgeindex1(refinededgeindex2), currincludededgeindex) = opvert;
                                    end
                                end
                                
                            end
                        else
                            %                     disp('EDGE NOT RELAXED')
                            %                     disp(edge);
                            %                     Sr(rcounter,:) = sort([edge(1) edge(2) size(vertices2,1)]);
                            %                     rcounter = rcounter + 1;
                            %                     % HAPON UPDATES
                            %                     opvert = setdiff(index, edge);
                            %                     curroppvertex = newinterioroppositevertex1(edgeindex1(edgeindex2,:),:);
                            %                     currincludededgeindex = find(curroppvertex == opvert);
                            %
                            %                     newinterioroppositevertex1(edgeindex1(edgeindex2,:),currincludededgeindex) = size(vertices2,1);
                            disp('EDGE NOT RELAXED')
                            disp(edge);
                            Sr(rcounter,:) = sort([edge(1) edge(2) size(vertices2,1)]);
                            rcounter = rcounter + 1;
                            
                            % HAPON UPDATES
                            opvert = setdiff(index, edge);
                            curroppvertex = newinterioroppositevertex1(edgeindex1(edgeindex2,:),:);
                            curroppvertex2 = newoppositevertex1(outervertindex1(outervertindex2,:),:);
                            
                            
                            currincludededgeindex = find(curroppvertex == opvert);
                            currincludededgeindex2 = find(curroppvertex2 == opvert);
                            
                            newinterioroppositevertex1(edgeindex1(edgeindex2,:),currincludededgeindex) = size(vertices2,1);
                            newoppositevertex1(outervertindex1(outervertindex2,:),currincludededgeindex2) = size(vertices2,1);
                        end
                        
                    end
                    
                    neighbors1 = find(newedges1(:,1) == size(vertices2,1));
                    neighbors2 = find(newedges1(:,2) == size(vertices2,1));
                    neighbors3 = find(newinterioredges1(:,1) == size(vertices2,1));
                    neighbors4 = find(newinterioredges1(:,2) == size(vertices2,1));
                    %disp([neighbors1; neighbors2; neighbors3; neighbors4]);
                    if (size([neighbors1; neighbors2; neighbors3; neighbors4],1) ~= 0),
                        aveedgelength(size(vertices2,1),:) = mean(newedgeslength([neighbors1; neighbors2; neighbors3; neighbors4],1));
                    end
                    break;
                else
                    if k == 3,
                        Sr(rcounter,:) = sort(index);
                        rcounter = rcounter + 1;
                    end
                end
            end
        end
        Sr = unique(Sr, 'rows');
        
        
        counter = 1;
        rowstodelete = [];
        for i = 1:size(Sr,1),
            row = Sr(i,:);
            if (length(unique(row)) ~= 3),
                rowstodelete(counter) = i;
                counter = counter + 1;
            end
        end
        
        Sr(rowstodelete,:) = [];
        
        patchedges = firsthole;                           %% Saves the patch edges per cluster
        for i = 1:size(patchedges,1),
            curredge = patchedges(i,:);
            x = line(vertices([curredge(1) curredge(2)],1),vertices([curredge(1) curredge(2)],2), vertices([curredge(1) curredge(2)],3));
            set(x, 'LineWidth', 5);
            set(x, 'Color', 'r');
        end
        
        patchedges = secondhole;                           %% Saves the patch edges per cluster
        for i = 1:size(patchedges,1),
            curredge = patchedges(i,:);
            x = line(vertices([curredge(1) curredge(2)],1),vertices([curredge(1) curredge(2)],2), vertices([curredge(1) curredge(2)],3));
            set(x, 'LineWidth', 5);
            set(x, 'Color', 'r');
        end
        
        
        
        
        
        
        
        
        
        
        
        finalfaces = [finalfaces; Sr];
        
%         for i = 1:size(Sr,1),
%             fprintf(fid, 'f ');
%             for j = 1:size(Sr,2),
%                 fprintf(fid, num2str(Sr(i,j)));
%                 if j<3,
%                     fprintf(fid, ' ');
%                 else
%                     fprintf(fid, '\n');
%                 end
%             end
%         end


    end
end



fid = fopen(['Output/CaptainAmerica.obj'], 'w');

for i = 1:size(vertices2,1),
    fprintf(fid, 'v ');
    for j = 1:size(vertices2,2),
        fprintf(fid, num2str(vertices2(i,j)));
        if j<3,
            fprintf(fid, ' ');
        else
            fprintf(fid, '\n');
        end
    end
end

Sr = [finalfaces; connectingFaces];
for i = 1:size(Sr,1),
    fprintf(fid, 'f ');
    for j = 1:size(Sr,2),
        fprintf(fid, num2str(Sr(i,j)));
        if j<3,
            fprintf(fid, ' ');
        else
            fprintf(fid, '\n');
        end
    end
end

fclose(fid);

toc

save(['ProcessedAll_', fName]);
save(['MeshCut_', fName], 'Sr', 'vertices2');
save(['MeshCutPaths_', fName], 'allPaths');