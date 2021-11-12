polytracks = [4];
fitpoly = zeros(n_tracks,1);
fitpoly(polytracks) = 1;
for j = 1:n_tracks
if fitpoly(j,1)>0
track_hwd{j} = [];
track_x{j} = [];
track_time{j} = [];
else
end
end
