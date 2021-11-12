%%Code to align kymograph base on the tracking of one protomer found using
%%S0_LS_Image analysis 
% After running this code re-run S0_LS_Image analysis 

track = 1; %align based on track 

%%
mx_shift = max(track_x{track})-min(track_x{track});
As = zeros(sz(1)+mx_shift+1,sz(2));

 shift =  max(track_x{track})-track_x{track} +1;

for i = 1:t_n
    As(shift(i):shift(i)+sz(1)-1,i) = A(1:sz(1),i);
end
kymo = 1;
 As = As-min(As(:));