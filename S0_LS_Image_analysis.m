%% Code to analyze kymograph images to output height vs time data for each high feature
%left to right should equal the time domain

%This code was developed in the Scheuring-Lab and used to analyze the motions 
%of GltPH transporters acquired by high-speed AFM line scanning. 
%% If using this code please refer to and cite:
%  Matin, T. R., Heath, G. R., Huysmans, G. H., Boudker, O., & Scheuring, S. 
%  (2020). Millisecond dynamics of an unlabeled amino acid transporter.
%  Nature communications, 11(1), 1-11. https://doi.org/10.1038/s41467-020-18811-z

%% Step 1
%run code - check output image that features have been tracked correctly,
%if not then alter below setting to optimize tracking. 
%For lateral drift compensation use S1_1_Kymo_align after running this code
%then re-run this code to analyze the drift corrected kymograph.
%If tracks are needed to be removed use S1_2_track_remover.

%% Step 2
%Go to S1_LS_Track_Analyzer 

%% detection and tracking variables:
avx = 20;           %pixel average in x (time domain)
avy = 7;             %pixel average in y (space domain)
thresh = 0.5;        %ignore traces with heights below thresh % of maximum

%tracking variables:
max_g = 2000;        %maximum down time (pix) (for tracking)
max_drift = 10;       %maxium +/- pixels to search for in x

%exclude edge or short tracks:
edg = 5;            %exclude protomers with edg pixels
min_length = 100;    %minimum protomers track length (pixels)


%% 
if exist('kymo')==0 
   if exist('f')==1
prompt = 'Open new image (0 for no, 1 for yes)? ';
new = input(prompt);
   else new=1;
   end
   
        if new ==1
        [f,path] = uigetfile('*.tif');
        f = fullfile(path,f); %tif (in nm) filename
        end
end

if exist('kymo')==0
    clearvars -except max_g avx avy min_length edg max_drift f thresh kymo time t_n 

    A = imread(f);
    t_n = numel(A(1,:));  
    x_n = numel(A(:,1));
    sz = size(A);
    time =linspace(1,t_n,t_n)';

    else
        clearvars -except A As avx avy f time t_n max_g min_length edg max_drift thresh kymo
        A = As;
        x_n = numel(A(:,1));
        sz = size(A);
        time =linspace(1,t_n,t_n)';
end

Amx = movmean(A,avy,1);
Amx = movmean(Amx,avx,2);

Amt = Amx.*(Amx>(thresh*max(Amx(:))));

Aall = reshape(Amt, [1, x_n*t_n]); %combines data into single colunm

   [h,locs,wd] = findpeaks(Aall,1,'MinPeakProminence',0.01,...
                                    'MinPeakDistance',2,...
                                    'WidthReference','halfheight');
                                    %'Annotate','extents');

peaks_n = numel(wd);
temptime = zeros(1,peaks_n);
%reassemble
for i = 1:t_n
   for j = 1:peaks_n
        if (locs(j) > (i-1)*x_n) && locs(j) < i*x_n
        temptime(:,j) = i;
        end
    end
end
locs = locs - x_n*(temptime-1);

for i = 1:t_n
pos = temptime>(i-0.1) & temptime<(i+1);
struc = locs(pos)';
struc_wd = wd(pos)';
CC{i} = struc;
CC_wd{i} = struc_wd;
end
CC = CC';
CC_wd = CC_wd';

%Track peaks 
tracks = simpletracker(CC,'Method','Hungarian','MaxLinkingDistance',max_drift,'MaxGapClosing',max_g);

n_tracks = numel(tracks)

tracked = zeros(t_n,n_tracks);
trackedN = tracked;  tracked_wd = tracked; trackedN_wd = tracked;

for i_track = 1:n_tracks
    pos = (tracks{i_track,1});
    for j = 1:t_n
        if pos(j)>0
        tracked(j,i_track) = CC{j,1}(pos(j));         trackedN(j,i_track) = CC{j,1}(pos(j));
        tracked_wd(j,i_track) = CC_wd{j,1}(pos(j));   trackedN_wd(j,i_track) = CC_wd{j,1}(pos(j));
        else
        tracked(j,i_track) = 0;
        trackedN(j,i_track) = NaN;
        tracked_wd(j,i_track) = 0;
        trackedN_wd(j,i_track) = NaN;
        end
    end
end
 
y = tracked;

 for j = 1:n_tracks      %creat tracks with start-finish points and 
        for i = 1:t_n
            if sum(y(1:i,j))> 0.9 && sum(y(i:numel(y(:,1)),j))> 0.9
            yt(i,j) = 1;
            else
            yt(i,j) = 0;
            end
        end   
 end
for j = 1:n_tracks
t1 = time(yt(:,j)>0.5);
 
track_time{j} = t1;
track_x{j} = trackedN(t1,j);
track_x{j} = fillmissing(track_x{j},'movmean',max_g,'EndValues','nearest');
track_x{j} = round(movmean(track_x{j},100));
        
track_wd{j} = trackedN_wd(t1,j);
track_wd{j} = fillmissing(track_wd{j},'nearest');
track_wd{j} = movmean(track_wd{j},100);
        
track_xu{j} = track_x{j}+ round(track_wd{j}/5);
track_xd{j} = track_x{j}- round(track_wd{j}/5);
   
 end

  for j = 1:n_tracks  
     for t_i = 1:numel(track_time{j}) 
        if track_x{j}< x_n- edg & track_x{j} > edg & numel(track_x{j}) >10 & sum(tracked(:,j) > 0.5) > min_length %exclude boundary tracks and tracks shorter than 300 pixels
        
        track_hi{j}(t_i) = A(track_x{j}(t_i),track_time{j}(t_i));  
        track_hwd{j}(t_i) = mean(A(track_xd{j}(t_i):track_xu{j}(t_i),track_time{j}(t_i)))';  
        
        else
         track_x{j} = [];
         track_time{j} = [];
         track_wd{j} = [];
        end
     end  
  end  

%% plot kymo
figure('Position',[50 50 1400 700])
tiledlayout(2,2, 'Padding', 'none', 'TileSpacing', 'none');
nexttile
imagesc(A)
title('input image')
colormap(jet)
nexttile
imagesc(Amx)
title('filtered')
nexttile
imagesc(Amt)
hold on
plot(temptime,locs,'.')
title('Threshold + detections')
nexttile
imagesc(A)
hold on
for j = 1:n_tracks  
    %if track_x{j}>0
text(min(track_time{1,j})+50,min(track_x{1,j}),num2str(j),'Color','c')
plot(track_time{1,j},track_x{1,j},'LineWidth',1.5)
    %end
end
title('Tracked')
colormap(gray)

%% clear extra data
clearvars CC CC_wd h pos yt y wd trackedN trackedN_wd tracked_wd time Aall t1 ti i i_track j peaks_n struc struc_wd 


%% functions for simpletracker
%Jean-Yves Tinevez (2021). simpletracker (https://github.com/tinevez/simpletracker), GitHub.
%BSD 3-Clause License
function [ tracks adjacency_tracks A ] = simpletracker(points, varargin)
% SIMPLETRACKER  a simple particle tracking algorithm that can deal with gaps

% INPUT SYNTAX
%
% tracks = SIMPLETRACKER(points) rebuilds the tracks generated by the
% particle whose coordinates are in |points|. |points| must be a cell
% array, with one cell per frame considered. Each cell then contains the
% coordinates of the particles found in that frame in the shape of a
% |n_points x n_dim| double array, where |n_points| is the number of points
% in that frame (that can vary a lot from one frame to another) and |n_dim|
% is the dimensionality of the problem (1 for 1D, 2 for 2D, 3 for 3D,
% etc...).
%
% tracks = SIMPLETRACKER(points, KEY, VALUE, ...)  allows to pass extra
% parameters to configure the tracking process settings. Accepted KEYS &
% VALUES are:
%
% 'MaxLinkingDistance' - a positive number, by default Inifity.
% Defines a maximal distance for particle linking. Two particles will not
% be linked (even if they are the remaining closest pair) if their distance
% is larger than this value. By default, it is infinite, not preventing nay
% linking.
% 
% 'MaxGapClosing' - a positive integer, by default 3
% Defines a maximal frame distance in gap-closing. Frames further way than
% this value will not be investigated for gap closing. By default, it has
% the value of 3.
%
% 
% VERSION HISTORY
%
% * v1.0 - November 2011 - Initial release.
% * v1.1 - May 2012 - Solve memory problems for large number of points.
%                   - Considerable speed improvement using properly the
%                   sparse matrices.
%                   - Use the key/value pair syntax to configure the
%                   function.
% * v1.3 - August 2012 - Fix a severe bug thanks to Dave Cade
%
% Jean-Yves Tinevez < jeanyves.tinevez@gmail.com> November 2011 - 2012

    %% Parse arguments
    
    p = inputParser;
    defaultDebug                = false;
    defaultMaxGapClosing        = 3;
    defaultMaxLinkingDistance   = Inf;
    defaultMethod               = 'Hungarian';
    expectedMethods = { defaultMethod, 'NearestNeighbor' };
    
    p.addParamValue('Debug', defaultDebug, @islogical);
    p.addParamValue('MaxGapClosing', defaultMaxGapClosing, @isnumeric);
    p.addParamValue('MaxLinkingDistance', defaultMaxLinkingDistance, @isnumeric);
    p.addParamValue('Method', defaultMethod,...
         @(x) any(validatestring(x, expectedMethods)));
    
    p.parse( varargin{:} );
    
    debug                   = p.Results.Debug;
    max_gap_closing         = p.Results.MaxGapClosing;
    max_linking_distance    = p.Results.MaxLinkingDistance;
    method                  = p.Results.Method;
    
    %% Frame to frame linking
    
    if debug
       fprintf('Frame to frame linking using %s method.\n', method);
    end
    
    n_slices = numel(points);
    
    current_slice_index = 0;
    row_indices = cell(n_slices, 1);
    column_indices = cell(n_slices, 1);
    unmatched_targets = cell(n_slices, 1);
    unmatched_sources = cell(n_slices, 1);
    n_cells = cellfun(@(x) size(x, 1), points);
    
    if debug
       fprintf('%03d/%03d', 0, n_slices-1);
    end
    
    for i = 1 : n_slices-1
        
        if debug
            fprintf(repmat('\b', 1, 7)); 
            fprintf('%03d/%03d', i, n_slices-1);
        end

        source = points{i};
        target = points{i+1};
        
        % Frame to frame linking
        switch lower(method)
        
            case 'hungarian'
                [target_indices , ~, unmatched_targets{i+1} ] = ...
                    hungarianlinker(source, target, max_linking_distance);
        
            case 'nearestneighbor'
                [target_indices , ~, unmatched_targets{i+1} ] = ...
                    nearestneighborlinker(source, target, max_linking_distance);

        end
        
        
        unmatched_sources{i} = find( target_indices == -1 );
        
        % Prepare holders for links in the sparse matrix
        n_links = sum( target_indices ~= -1 );
        row_indices{i} = NaN(n_links, 1);
        column_indices{i} = NaN(n_links, 1);
        
        % Put it in the adjacency matrix
        index = 1;
        for j = 1 : numel(target_indices)
            
            % If we did not find a proper target to link, we skip
            if target_indices(j) == -1
                continue
            end
            
            % The source line number in the adjacency matrix
            row_indices{i}(index) = current_slice_index + j;
            
            % The target column number in the adjacency matrix
            column_indices{i}(index) = current_slice_index + n_cells(i) + target_indices(j);
            
            index = index + 1;
            
        end
        
        current_slice_index = current_slice_index + n_cells(i);
        
    end

 
    row_index = vertcat(row_indices{:});
    column_index = vertcat(column_indices{:});
    link_flag = ones( numel(row_index), 1);
    n_total_cells = sum(n_cells);
    
    if debug
        fprintf('\nCreating %d links over a total of %d points.\n', numel(link_flag), n_total_cells)
    end

    A = sparse(row_index, column_index, link_flag, n_total_cells, n_total_cells);
    
    if debug
        fprintf('Done.\n')
    end
    
    
    %% Gap closing
    
    if debug
        fprintf('Gap-closing:\n')
    end
    
    current_slice_index = 0;
    for i = 1 : n_slices-2
        
        
        % Try to find a target in the frames following, starting at i+2, and
        % parsing over the target that are not part in a link already.
        
        current_target_slice_index = current_slice_index + n_cells(i) + n_cells(i+1);
        
        for j = i + 2 : min(i +  max_gap_closing, n_slices)
            
            source = points{i}(unmatched_sources{i}, :);
            target = points{j}(unmatched_targets{j}, :);
            
            if isempty(source) || isempty(target)
                current_target_slice_index = current_target_slice_index + n_cells(j);
                continue
            end
            
            target_indices = nearestneighborlinker(source, target, max_linking_distance);
            
            % Put it in the adjacency matrix
            for k = 1 : numel(target_indices)
                
                % If we did not find a proper target to link, we skip
                if target_indices(k) == -1
                    continue
                end
                
                if debug
                    fprintf('Creating a link between point %d of frame %d and point %d of frame %d.\n', ...
                        unmatched_sources{i}(k), i, unmatched_targets{j}(target_indices(k)), j);
                end
                
                % The source line number in the adjacency matrix
                row_index = current_slice_index + unmatched_sources{i}(k);
                % The target column number in the adjacency matrix
                column_index = current_target_slice_index + unmatched_targets{j}(target_indices(k));
                
                A(row_index, column_index) = 1; %#ok<SPRIX>
                
            end
            
            new_links_target =  target_indices ~= -1 ;
            
            % Make linked sources unavailable for further linking
            unmatched_sources{i}( new_links_target ) = [];
            
            % Make linked targets unavailable for further linking
            unmatched_targets{j}(target_indices(new_links_target)) = [];
            
            current_target_slice_index = current_target_slice_index + n_cells(j);
        end
        
        current_slice_index = current_slice_index + n_cells(i);
        
    end
    
    if debug
        fprintf('Done.\n')
    end
    
    %% Parse adjacency matrix to build tracks
    
    if debug
        fprintf('Building tracks:\n')
    end
    
    % Find columns full of 0s -> means this cell has no source
    cells_without_source = [];
    for i = 1 : size(A, 2)
        if length(find(A(:,i))) == 0 %#ok<ISMT>
            cells_without_source = [ cells_without_source ; i ]; %#ok<AGROW>
        end
    end
    
    n_tracks = numel(cells_without_source);
    adjacency_tracks = cell(n_tracks, 1);
    
    AT = A';
    
    for i = 1 : n_tracks
        
        tmp_holder = NaN(n_total_cells, 1);
        
        target = cells_without_source(i);
        index = 1;
        while ~isempty(target)
            tmp_holder(index) = target;
            target = find( AT(:, target), 1, 'first' );
            index = index + 1;
        end
        
        adjacency_tracks{i} = tmp_holder ( ~isnan(tmp_holder) );
    end
    
    %% Reparse adjacency track index to have it right.
    % The trouble with the previous track index is that the index in each
    % track refers to the index in the adjacency matrix, not the point in
    % the original array. We have to reparse it to put it right.
    
    tracks = cell(n_tracks, 1);
    
    for i = 1 : n_tracks
        
        adjacency_track = adjacency_tracks{i};
        track = NaN(n_slices, 1);
        
        for j = 1 : numel(adjacency_track)
            
            cell_index = adjacency_track(j);
            
            % We must determine the frame this index belong to
            tmp = cell_index;
            frame_index = 1;
            while tmp > 0
                tmp = tmp - n_cells(frame_index);
                frame_index = frame_index + 1;
            end
            frame_index = frame_index - 1;
            in_frame_cell_index = tmp + n_cells(frame_index);
            
            track(frame_index) = in_frame_cell_index;
            
        end
        
        tracks{i} = track;
        
    end
    
end
function [ target_indices target_distances unassigned_targets total_cost ] = hungarianlinker(source, target, max_distance)
%HUNGARIANLINKER link two lists of points based on the hungarian algorithm.
%
% target_indices = HUNGARIANLINKER(source, target) finds for each point in
% 'source' the closest point in 'target'. These 2 inputs must be arrays
% with one point per row, and have their cartesian coordinates in each
% column (1D, 2D, 3D, ...). Source to target assignment is based on the
% famous hungarian algorithm using its excellent implementation by the
% excellent Yi Cao. The two arrays might not have the same number of
% points.
%
% The indices of the 'target' points are returned in an array
% 'target_indices', so that each row in 'source' matches the corresponding
% row in 'target(target_indices, :)'.
%
% The linking is exclusive: one source point is linked to at most one
% target point, and conversely. The linking is globally optimal: the sum of
% the square distance is minimized, contrary to the naive nearest neighbor
% approach.
%
% target_indices = HUNGARIANLINKER(source, target, max_distance) adds
% a condition on distance. If the nearest neighbor is found to be at a
% distance larger than the given 'max_distance', they are not linked, and
% the 'target_indices' receive the value -1 for this source point. The same
% happens if all target points are exhausted.
% 
% [ target_indices target_distances ] = HUNGARIANLINKER(source, target)
% additionaly return the distance to the matched target point. Un-matched
% source points have a distance value set to NaN.
%
% [ target_indices target_distances unmatched_targets ] =
%                                         HUNGARIANLINKER(source, target) 
% additionaly return the indices of the points in 'target' that have not
% been linked.
%
% [ target_indices target_distances unmatched_targets total_cost ] =
%                                         HUNGARIANLINKER(source, target) 
% additionaly return the globally optimized value of the square distance
% sum.
%
% The matching algorithm used here is one of the best available and ensures
% that the resulting assignment is a optimum. However the price to pay is
% an increased complexity. The cost for the naive nearest neighbor approach
% roughly scales as O(p^2) where p is the number of source points. The
% munkres implementation of the hungarian algorithm by Yi Cao is in O(p^3),
% and is the best so far.
%
% EXAMPLE:
% 
% n_points = 20;
% source = 10 * rand(n_points, 2);
% target = source + rand(n_points, 2);
% target_indices = hungarianlinker(source, target);
% colors = hsv(n_points);
% figure
% hold on
% for i = 1 :n_points
%    plot(source(i,1), source(i,2), 'o', 'Color', colors(i,:))
%    plot(target(target_indices(i),1), target(target_indices(i),2), 's', ...
%       'Color', colors(i,:))
%    plot( [ source(i,1) target(target_indices(i),1) ] , ...
%       [ source(i,2)  target(target_indices(i),2) ], ...
%        'Color', colors(i,:))
% end
% 
%
% Jean-Yves Tinevez <jeanyves.tinevez@gmail.com>.
% However all credits should go to Yi Cao, which did the hard job of
% implementing the Munkres algorithm; this file is merely a wrapper for it.

    if nargin < 3
        max_distance = Inf;
    end

    n_source_points = size(source, 1);
    n_target_points = size(target, 1);
    
    D = NaN(n_source_points, n_target_points);
    
    % Build distance matrix
    for i = 1 : n_source_points
        
        % Pick one source point
        current_point = source(i, :);
        
        % Compute square distance to all target points
        diff_coords = target - repmat(current_point, n_target_points, 1);
        square_dist = sum(diff_coords.^2, 2);
        
        % Store them
        D(i, :) = square_dist;
        
    end
    
    % Deal with maximal linking distance: we simply mark these links as already
    % treated, so that they can never generate a link.
    D ( D > max_distance * max_distance ) = Inf;
    
    % Find the optimal assignment is simple as calling Yi Cao excellent FEX
    % submission.
    [ target_indices total_cost ] = munkres(D);
    % Set unmatched sources to -1
    target_indices ( target_indices  == 0 ) = -1;
    
    % Collect distances
    target_distances = NaN(numel(target_indices), 1);
    for i = 1 : numel(target_indices)
        if target_indices(i) < 0
            continue
        end
        
        target_distances(i) = sqrt ( D ( i , target_indices(i)) );
        
    end
    
    unassigned_targets = setdiff ( 1 : n_target_points , target_indices );
    
    
end
function [assignment,cost] = munkres(costMat)
% MUNKRES   Munkres (Hungarian) Algorithm for Linear Assignment Problem. 
%
% [ASSIGN,COST] = munkres(COSTMAT) returns the optimal column indices,
% ASSIGN assigned to each row and the minimum COST based on the assignment
% problem represented by the COSTMAT, where the (i,j)th element represents the cost to assign the jth
% job to the ith worker.
%
% Partial assignment: This code can identify a partial assignment is a full
% assignment is not feasible. For a partial assignment, there are some
% zero elements in the returning assignment vector, which indicate
% un-assigned tasks. The cost returned only contains the cost of partially
% assigned tasks.

% This is vectorized implementation of the algorithm. It is the fastest
% among all Matlab implementations of the algorithm.

% Examples
% Example 1: a 5 x 5 example
%{
[assignment,cost] = munkres(magic(5));
disp(assignment); % 3 2 1 5 4
disp(cost); %15
%}
% Example 2: 400 x 400 random data
%{
n=400;
A=rand(n);
tic
[a,b]=munkres(A);
toc                 % about 2 seconds 
%}
% Example 3: rectangular assignment with inf costs
%{
A=rand(10,7);
A(A>0.7)=Inf;
[a,b]=munkres(A);
%}
% Example 4: an example of partial assignment
%{
A = [1 3 Inf; Inf Inf 5; Inf Inf 0.5]; 
[a,b]=munkres(A)
%}
% a = [1 0 3]
% b = 1.5
% Reference:
% "Munkres' Assignment Algorithm, Modified for Rectangular Matrices", 
% http://csclab.murraystate.edu/bob.pilgrim/445/munkres.html

% version 2.3 by Yi Cao at Cranfield University on 11th September 2011

assignment = zeros(1,size(costMat,1));
cost = 0;

validMat = costMat == costMat & costMat < Inf;
bigM = 10^(ceil(log10(sum(costMat(validMat))))+1);
costMat(~validMat) = bigM;

% costMat(costMat~=costMat)=Inf;
% validMat = costMat<Inf;
validCol = any(validMat,1);
validRow = any(validMat,2);

nRows = sum(validRow);
nCols = sum(validCol);
n = max(nRows,nCols);
if ~n
    return
end

maxv=10*max(costMat(validMat));

dMat = zeros(n) + maxv;
dMat(1:nRows,1:nCols) = costMat(validRow,validCol);

%*************************************************
% Munkres' Assignment Algorithm starts here
%*************************************************

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   STEP 1: Subtract the row minimum from each row.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
minR = min(dMat,[],2);
minC = min(bsxfun(@minus, dMat, minR));

%**************************************************************************  
%   STEP 2: Find a zero of dMat. If there are no starred zeros in its
%           column or row start the zero. Repeat for each zero
%**************************************************************************
zP = dMat == bsxfun(@plus, minC, minR);

starZ = zeros(n,1);
while any(zP(:))
    [r,c]=find(zP,1);
    starZ(r)=c;
    zP(r,:)=false;
    zP(:,c)=false;
end

while 1
%**************************************************************************
%   STEP 3: Cover each column with a starred zero. If all the columns are
%           covered then the matching is maximum
%**************************************************************************
    if all(starZ>0)
        break
    end
    coverColumn = false(1,n);
    coverColumn(starZ(starZ>0))=true;
    coverRow = false(n,1);
    primeZ = zeros(n,1);
    [rIdx, cIdx] = find(dMat(~coverRow,~coverColumn)==bsxfun(@plus,minR(~coverRow),minC(~coverColumn)));
    while 1
        %**************************************************************************
        %   STEP 4: Find a noncovered zero and prime it.  If there is no starred
        %           zero in the row containing this primed zero, Go to Step 5.  
        %           Otherwise, cover this row and uncover the column containing 
        %           the starred zero. Continue in this manner until there are no 
        %           uncovered zeros left. Save the smallest uncovered value and 
        %           Go to Step 6.
        %**************************************************************************
        cR = find(~coverRow);
        cC = find(~coverColumn);
        rIdx = cR(rIdx);
        cIdx = cC(cIdx);
        Step = 6;
        while ~isempty(cIdx)
            uZr = rIdx(1);
            uZc = cIdx(1);
            primeZ(uZr) = uZc;
            stz = starZ(uZr);
            if ~stz
                Step = 5;
                break;
            end
            coverRow(uZr) = true;
            coverColumn(stz) = false;
            z = rIdx==uZr;
            rIdx(z) = [];
            cIdx(z) = [];
            cR = find(~coverRow);
            z = dMat(~coverRow,stz) == minR(~coverRow) + minC(stz);
            rIdx = [rIdx(:);cR(z)];
            cIdx = [cIdx(:);stz(ones(sum(z),1))];
        end
        if Step == 6
            % *************************************************************************
            % STEP 6: Add the minimum uncovered value to every element of each covered
            %         row, and subtract it from every element of each uncovered column.
            %         Return to Step 4 without altering any stars, primes, or covered lines.
            %**************************************************************************
            [minval,rIdx,cIdx]=outerplus(dMat(~coverRow,~coverColumn),minR(~coverRow),minC(~coverColumn));            
            minC(~coverColumn) = minC(~coverColumn) + minval;
            minR(coverRow) = minR(coverRow) - minval;
        else
            break
        end
    end
    %**************************************************************************
    % STEP 5:
    %  Construct a series of alternating primed and starred zeros as
    %  follows:
    %  Let Z0 represent the uncovered primed zero found in Step 4.
    %  Let Z1 denote the starred zero in the column of Z0 (if any).
    %  Let Z2 denote the primed zero in the row of Z1 (there will always
    %  be one).  Continue until the series terminates at a primed zero
    %  that has no starred zero in its column.  Unstar each starred
    %  zero of the series, star each primed zero of the series, erase
    %  all primes and uncover every line in the matrix.  Return to Step 3.
    %**************************************************************************
    rowZ1 = find(starZ==uZc);
    starZ(uZr)=uZc;
    while rowZ1>0
        starZ(rowZ1)=0;
        uZc = primeZ(rowZ1);
        uZr = rowZ1;
        rowZ1 = find(starZ==uZc);
        starZ(uZr)=uZc;
    end
end

% Cost of assignment
rowIdx = find(validRow);
colIdx = find(validCol);
starZ = starZ(1:nRows);
vIdx = starZ <= nCols;
assignment(rowIdx(vIdx)) = colIdx(starZ(vIdx));
pass = assignment(assignment>0);
pass(~diag(validMat(assignment>0,pass))) = 0;
assignment(assignment>0) = pass;
cost = trace(costMat(assignment>0,assignment(assignment>0)));

function [minval,rIdx,cIdx]=outerplus(M,x,y)
ny=size(M,2);
minval=inf;
for c=1:ny
    M(:,c)=M(:,c)-(x+y(c));
    minval = min(minval,min(M(:,c)));
end
[rIdx,cIdx]=find(M==minval);
end
end
function [ target_indices target_distances unassigned_targets ] = nearestneighborlinker(source, target, max_distance)
%NEARESTNEIGHBORLINKER link two lists of points based on nearest neighbor.
%
% target_indices = NEARESTNEIGHBORLINKER(source, target) finds for each
% point in 'source' the closest point in 'target'. These 2 inputs must be
% arrays with one point per row, and have their cartesian coordinates in
% each column (1D, 2D, 3D, ...). Nearest neighbor matching is based on
% euclidean distance. The two arrays might not have the same number of
% points.
%
% The indices of the 'target' points are returned in an array
% 'target_indices', so that each row in 'source' matches the corresponding
% row in 'target(target_indices, :)'.
%
% The linking is exclusive: one source point is linked to at most one
% target point, and conversely. The linking is only locally optimal: the
% two closest points amongst the two sets are sought for first, then the
% second closest pair, excluding the first, etc... This ensures that the
% resulting linking will not depend on the order of the points in each set.
%
% target_indices = NEARESTNEIGHBORLINKER(source, target, max_distance) adds
% a condition on distance. If the nearest neighbor is found to be at a
% distance larger than the given 'max_distance', they are not linked, and
% the 'target_indices' receive the value -1 for this source point. The same
% happens if all target points are exhausted.
% 
% [ target_indices target_distances ] = 
%                                   NEARESTNEIGHBORLINKER(source, target)
% additionaly return the distance to the matched target point. Un-matched
% source points have a distance value set to NaN.
%
% [ target_indices target_distances unmatched_targets ]= 
%                                   NEARESTNEIGHBORLINKER(source, target)
% additionaly return the indices of the points in 'target' that have not
% been linked.
%
% This is the cheapest (in term of accuracy) algorithm for linking that can
% be made. In particular, it is not guaranteed (and it is generally not the
% case) that the returned linking is an optimum for the sum of distances.
% Each source point is matched regardless of the others, there is no global
% optimization here (the Hungarian algorithm does that). Also, there exists
% refinement to nearest neighbor searches, such as the use of KD-trees;
% this contribution is exempt of such developments.
%
% EXAMPLE:
% 
% n_points = 20;
% source = 10 * rand(n_points, 2);
% target = source + rand(n_points, 2);
% target_indices = nearestneighborlinker(source, target);
% colors = hsv(n_points);
% figure
% hold on
% for i = 1 :n_points
%    plot(source(i,1), source(i,2), 'o', 'Color', colors(i,:))
%    plot(target(target_indices(i),1), target(target_indices(i),2), 's', ...
%       'Color', colors(i,:))
%    plot( [ source(i,1) target(target_indices(i),1) ] , ...
%       [ source(i,2)  target(target_indices(i),2) ], ...
%        'Color', colors(i,:))
% end
% 
% VERSION HISTORY
%
% * v1.0 - November 2011 - Initial release.
% * v1.1 - May 2012 - Fix a severe bug thanks to Dave Cade
%
% Jean-Yves Tinevez < jeanyves.tinevez@gmail.com> November 2011 - 2012

    if nargin < 3
        max_distance = Inf;
    end
   
    n_source_points = size(source, 1);
    n_target_points = size(target, 1);
    
    D = NaN(n_source_points, n_target_points);
    
    % Build distance matrix
    for i = 1 : n_source_points
        
        % Pick one source point
        current_point = source(i, :);
        
        % Compute square distance to all target points
        diff_coords = target - repmat(current_point, n_target_points, 1);
        square_dist = sum(diff_coords.^2, 2);
        
        % Store them
        D(i, :) = square_dist;
        
    end
    
    % Deal with maximal linking distance: we simply mark these links as already
    % treated, so that they can never generate a link.
    D ( D > max_distance * max_distance ) = Inf;
    
    target_indices = -1 * ones(n_source_points, 1);
    target_distances = NaN(n_source_points, 1);
    
    % Parse distance matrix
    while ~all(isinf(D(:)))
        
        [ min_D closest_targets ] = min(D, [], 2); % index of the closest target for each source points
        [ ~, sorted_index ] = sort(min_D);
        
        for i = 1 : numel(sorted_index)
            
            source_index =  sorted_index(i);
            target_index =  closest_targets ( sorted_index(i) );
            
            % Did we already assigned this target to a source?
            if any ( target_index == target_indices )
                
                % Yes, then exit the loop and change the distance matrix to
                % prevent this assignment
                break
                
            else
                
                % No, then store this assignment
                target_indices( source_index ) = target_index;
                target_distances ( source_index ) = sqrt ( min_D (  sorted_index(i) ) );
                
                % And make it impossible to find it again by putting the target
                % point to infinity in the distance matrix
                D(:, target_index) = Inf;
                % And the same for the source line
                D(source_index, :) = Inf;
                
                if all(isinf(D(:)))
                    break
                end
                
            end
            
        end
        
    end
    
    unassigned_targets = setdiff ( 1 : n_target_points , target_indices );
        
end