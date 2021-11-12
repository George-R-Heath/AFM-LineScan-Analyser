%% Code to perform state analysis of height vs time traces extracted from kymograph
%input polynomial order into poly_all to apply a polynomial fitting to z
%data to all tracks. Use poly and polytracks to input polynomial to apply
%polynomial fitting to a single track.
%Use 'Quality' setting to exclude tracks with low signal to noise ratio

%dwell times for all traces output as 'all_down_wd', 'all_mid_wd' and 'all_up_wd'
%dwell heights for all traces output as 'all_down_h', 'all_mid_h' and 'all_up_h'
%normalized dwell heights for all traces output as 'all_down_h_n', 'all_mid_h_n' and 'all_up_h_n'

%           Important Output Values:
%  1. Track stats table for values for each track
%  2. s2 and s3 individual 2 and 3 state fits (depedning on state number given in
%     track stats) as well as dwell times in each state 
%  3. all_down_h, all_mid_h, all_up_h dwell heights in down, mid and up states
%  4. all_down_h_n... normalized dwell height in down, mid and up states
%  5. all_down_wd... dwell times in down, mid and up states

%This codes uses the STaSI method (please cite DOI: 10.1021/jz501435p)
%The original STaSI code included here is also avaliable at: https://github.com/LandesLab/STaSI

%% variables to modify:
Quality = 1.35;    %Quality threshold -> use only tracks with a (upper/lower state difference)/(track standard deviation) > Quality
stasi = 1;         %perform stasi?

poly_all = 0;      %perform poly_all oreder polynomial fit to all tracks
poly = 0;          %selective polynomial fit to certain tracks
polytracks = [];   %select tracks for poly fit

%% analysis code:
clearvars Amx Amt locs temptime

fitpoly = zeros(n_tracks,1);
fitpoly(polytracks) = 1;

if stasi == 1
    for j = 1:n_tracks 
      if track_x{j}< x_n-edg & sum(track_x{j} >edg)>0 & numel(track_x{j}) >5 & sum(tracked(:,j) > 0.5) > 300
      
          for i = 1:numel(track_hwd{j})-1
          delt_h{j}(i) = track_hwd{j}(i)-track_hwd{j}(i+1);
          end
      track_sd{j} = std(delt_h{j}); 
       
      [p,sii,mu] = polyfit((1:numel(track_hwd{j})'),track_hwd{j},poly_all);     %polynomial fit
      f_z = polyval(p,(1:numel(track_hwd{j}))',[],mu);
      track_hwd_p{j} = (track_hwd{j}' -f_z+mean(f_z))'; 
      output(j,:) = StaSI_inputeff(track_hwd_p{j}, track_sd{j});%polynomial fit subtractmin_h(j) = min(track_hwd{j});
      %output(j,:) = StaSI_inputeff(track_hwd{j}, track_sd{j});   
      max_h(j,1) = max(track_hwd_p{j});
      var_h(j,1) = var(track_hwd_p{j});
    
      end
    end
err_mult =1.4;   
for j = 1:n_tracks
if fitpoly(j,1)>0
      [p,sii,mu] = polyfit((1:numel(track_hwd{j})'),track_hwd{j},poly);     %polynomial fit
      f_z = polyval(p,(1:numel(track_hwd{j}))',[],mu);
      track_hwd_p{j} = (track_hwd{j}' -f_z +mean(f_z))'; 
      output(j,:) = StaSI_inputeff(track_hwd_p{j}, track_sd{j});%polynomial fit subtractmin_h(j) = min(track_hwd{j});
      %output(j,:) = StaSI_inputeff(track_hwd{j}, track_sd{j});   
      max_h(j,1) = max(track_hwd_p{j});
      var_h(j,1) = var(track_hwd_p{j});
      
else
end
end 
    
    
    for i = 1:numel(output)
        if output(i).breaks >0.5
            
            [minMDL(i,2),minMDL(i,1)] = min(output(i).MDL);
        end
    end
end
s2 = [];
s3 = [];
for j = 1:numel(output) 
  if output(j).breaks >0.5 & numel(track_x{j}) & numel(track_hwd_p{j}) 
    %2state   
stasi_2{j} = output(j).eff_fit(2,:);
[s2{j}.up_h, s2{j}.up_wd, s2{j}.down_h, s2{j}.down_wd, s2{j}.ideal, s2{j}.ideal_i]=TwoState_measure(track_time{j}, track_hwd_p{j},stasi_2{j});
[s2{j}.up_h, s2{j}.up_wd, s2{j}.down_h, s2{j}.down_wd, s2{j}.ideal, s2{j}.ideal_i]=TwoState_measure(track_time{j}, track_hwd_p{j},s2{j}.ideal);    
[s2{j}.up_h, s2{j}.up_wd, s2{j}.down_h, s2{j}.down_wd, s2{j}.ideal, s2{j}.ideal_i]=TwoState_measure(track_time{j}, track_hwd_p{j},s2{j}.ideal);  

%3state  
stasi_3{j} = output(j).eff_fit(3,:);
[s3{j}.up_h, s3{j}.up_wd, s3{j}.mid_h, s3{j}.mid_wd, s3{j}.down_h, s3{j}.down_wd, s3{j}.ideal, s3{j}.ideal_i]=ThreeState_measure(track_time{j}, track_hwd_p{j},stasi_3{j});
[s3{j}.up_h, s3{j}.up_wd, s3{j}.mid_h, s3{j}.mid_wd, s3{j}.down_h, s3{j}.down_wd, s3{j}.ideal, s3{j}.ideal_i]=ThreeState_measure(track_time{j}, track_hwd_p{j},s3{j}.ideal);
[s3{j}.up_h, s3{j}.up_wd, s3{j}.mid_h, s3{j}.mid_wd, s3{j}.down_h, s3{j}.down_wd, s3{j}.ideal, s3{j}.ideal_i]=ThreeState_measure(track_time{j}, track_hwd_p{j},s3{j}.ideal);

[states(j,1), states(j,2)] = State_numbers(track_hwd_p{j}, s2{j}.ideal, s3{j}.ideal,err_mult);
states(j,3) = (max(output(j).eff_fit(2,:)))  -  (min(output(j).eff_fit(2,:)));
states(j,4) = var(track_hwd_p{j});
[states(j,6) ,states(j,5)] = min(output(j).MDL);
states(j,7) = j;
states(j,8) = track_sd{j};
states(j,9) = states(j,3)/track_sd{j};
  end
end

s2 = [];
s3 = [];
for j = 1:numel(output)
  if states(j,2)>0 && output(j).breaks >0.5

    %2state   
stasi_2{j} = output(j).eff_fit(2,:);
stasi_3{j} = output(j).eff_fit(3,:);

states(j,3) = (max(output(j).eff_fit(2,:)))  -  (min(output(j).eff_fit(2,:)));         
for i = 1:numel(track_hwd{j})-1
          delt_h{j}(i) = track_hwd{j}(i)-track_hwd{j}(i+1);
end
track_sd{j} = std(delt_h{j}); 
states(j,8) = track_sd{j};
states(j,9) = states(j,3)/track_sd{j};

stasi_2_shift{j} = ThreeState_ToTwo(track_time{j}, track_hwd{j},stasi_3{j}); 
[s2{j}.up_h, s2{j}.up_wd, s2{j}.down_h, s2{j}.down_wd, s2{j}.ideal, s2{j}.ideal_i]=TwoState_measure(track_time{j}, track_hwd_p{j},stasi_2_shift{j});
[s2{j}.up_h, s2{j}.up_wd, s2{j}.down_h, s2{j}.down_wd, s2{j}.ideal, s2{j}.ideal_i]=TwoState_measure(track_time{j}, track_hwd_p{j},s2{j}.ideal);   
s2{j}.up_h =[]; s2{j}.up_wd =[];s2{j}.down_h =[]; s2{j}.down_wd =[];
[s2{j}.up_h, s2{j}.up_wd, s2{j}.down_h, s2{j}.down_wd, s2{j}.ideal, s2{j}.ideal_i]=TwoState_measure(track_time{j}, track_hwd_p{j},s2{j}.ideal);  
states(j,10) = mean((track_hwd_p{j}-s2{j}.ideal).^2).^0.5;
%3state  

[s3{j}.up_h, s3{j}.up_wd, s3{j}.mid_h, s3{j}.mid_wd, s3{j}.down_h, s3{j}.down_wd, s3{j}.ideal, s3{j}.ideal_i]=ThreeState_measure(track_time{j}, track_hwd_p{j},stasi_3{j});
[s3{j}.up_h, s3{j}.up_wd, s3{j}.mid_h, s3{j}.mid_wd, s3{j}.down_h, s3{j}.down_wd, s3{j}.ideal, s3{j}.ideal_i]=ThreeState_measure(track_time{j}, track_hwd_p{j},s3{j}.ideal);
[s3{j}.up_h, s3{j}.up_wd, s3{j}.mid_h, s3{j}.mid_wd, s3{j}.down_h, s3{j}.down_wd, s3{j}.ideal, s3{j}.ideal_i]=ThreeState_measure(track_time{j}, track_hwd_p{j},s3{j}.ideal);
 
  end
end

all_up_h  = [];all_down_h =[];all_up_wd  = [];all_down_wd =[];
all_mid_h = []; all_mid_wd = []; all_up_h_n = []; all_down_h_n =[]; all_mid_h_n = [];


for i = 1:numel(states(:,1))
%% 2 state analyis
if states(i,9) > Quality && states(i,2) ==2

            all_up_h = vertcat(all_up_h,s2{i}.up_h(:));
            all_up_h_n = vertcat(all_up_h_n,(s2{i}.up_h(:)-mean(s2{i}.down_h(:)))/states(i,3));
            all_up_wd = vertcat(all_up_wd,s2{i}.up_wd(:));
            
            all_down_h = vertcat(all_down_h,s2{i}.down_h(:));
            all_down_h_n = vertcat(all_down_h_n,(s2{i}.down_h(:)-mean(s2{i}.down_h(:)))/states(i,3));
            all_down_wd = vertcat(all_down_wd,s2{i}.down_wd(:));
     
end

        
%% 3 state analysis
 if states(i,9) > Quality && states(i,2) ==3
            all_up_h = vertcat(all_up_h,s3{i}.up_h(:));
            all_up_h_n = vertcat(all_up_h_n,(s3{i}.up_h(:)-mean(s3{i}.down_h(:)))/states(i,3));
            all_up_wd = vertcat(all_up_wd,s3{i}.up_wd(:));
            
            all_mid_h = vertcat(all_mid_h,s3{i}.mid_h(:));
            all_mid_h_n = vertcat(all_mid_h_n,(s3{i}.mid_h(:)-mean(s3{i}.down_h(:)))/states(i,3));
            all_mid_wd = vertcat(all_mid_wd,s3{i}.mid_wd(:));
            
            all_down_h = vertcat(all_down_h,s3{i}.down_h(:));
            all_down_h_n = vertcat(all_down_h_n,(s3{i}.down_h(:)-mean(s3{i}.down_h(:)))/states(i,3));
            all_down_wd = vertcat(all_down_wd,s3{i}.down_wd(:));
            
  end  
end
%% Figure plotting      

non_act = sum(states(:,2) == 1);
act = sum(states(:,2) > 1);

figure(7)
set(gcf, 'Name', 'active/non active','Position',[0 500 700 400])
imagesc(A)
colormap(gray)
hold on
for j = 1:numel(output) 
    if  states(j,2) >1
plot(track_time{1,j},track_x{1,j},'c','LineWidth',1)
text(min(track_time{1,j})+30,min(track_x{1,j}),num2str(j),'Color','c')
    else
text(min(track_time{1,j})+30,min(track_x{1,j}),num2str(j),'Color','r')
plot(track_time{1,j},track_x{1,j},'r')
    end
end
set(gca,'TickDir','out')


if act>0
figure(3) 
set(gcf, 'Name', '2 or 3 states','Position',[700 500 700 400])

ha = tight_subplot(2*act,1,[.0 .03],[.1 .01],[.1 .01]);
count = 0;
for track = 1:numel(output)
   if states(track,2) >1 & track_x{track}< x_n-edg & sum(track_x{track} >edg)>0
   count = count +1;   
    axes(ha(count));
    imagesc(A((track_x{track}-10):(track_x{track}+10),track_time{track}),[min(track_hwd_p{track}) max(track_hwd_p{track})]);
 
    set(gca,'TickDir','out')
    ylabel(num2str(track))
      hold on
   imagesc(0,20,1*(s2{track}.ideal>min(s2{track}.ideal)))
   set(gca,'TickDir','out')
   
    count = count +1;   
    axes(ha(count));
    imagesc(1*(s2{track}.ideal>min(s2{track}.ideal)))
    set(gca,'TickDir','out')
    hold on
    plot(track_hwd_p{track},'r','LineWidth',1)
    set(gca,'TickDir','out')
    if states(track,2) ==2
     plot(s2{track}.ideal,'b','LineWidth',1)
    else if states(track,2)==3
     plot(s3{track}.ideal,'g','LineWidth',2)
        end
    end
     xlim([0 numel(track_time{track})])
    ylim([min(track_hwd_p{track}) max(track_hwd_p{track})])
    set(ha(1),'XTickLabel','')
    set(gca,'Ydir','Normal')
    set(gca,'TickDir','out')
    end
end 
set(ha(1:(act-1)),'XTickLabel','')
end
colormap(gray)
if non_act>0
figure(4)
set(gcf, 'Name', '1 state','Position',[700 000 700 400])
ha = tight_subplot(2*non_act,1,[.0 .03],[.1 .01],[.1 .01]);
count = 0;
for track = 1:numel(output)
   if states(track,2) < 2 & track_x{track}< x_n-edg & sum(track_x{track} >edg)>0
      count = count +1;   
    axes(ha(count));
    imagesc(A((track_x{track}-10):(track_x{track}+10),track_time{track}),[min(track_hwd_p{track}) max(track_hwd_p{track})]);
   set(gca,'TickDir','out')
    ylabel(num2str(track))
      hold on
   imagesc(0,20,1*(s2{track}.ideal>min(s2{track}.ideal)))
 
    count = count +1;   
    axes(ha(count));
    imagesc(1*(s2{track}.ideal>min(s2{track}.ideal)))
    hold on
    plot(track_hwd_p{track})
        plot(output(track).eff_fit(1,:))
        set(gca,'TickDir','out')
    xlim([0 numel(track_time{track})])
    ylim([min(track_hwd_p{track}) max(track_hwd_p{track})])
    set(ha(1),'XTickLabel','')
    set(gca,'Ydir','Normal')
    end
end 
colormap(gray)
set(ha(1:(non_act-1)),'XTickLabel','')
end

figure(5)
set(gcf, 'Name', 'Dwell times','Position',[00 000 700 400])
plot(all_down_wd,all_down_h_n,'o')
hold on
plot(all_up_wd,all_up_h_n,'o')
plot(all_mid_wd,all_mid_h_n,'o')

set(gca,'XScale','Log','TickDir','out')
xlabel('Dwell time')
ylabel('Normalized dwell height')


Track_stats = array2table(states);
Track_stats.Properties.VariableNames = {'states' 'states (noise adjusted)' 'max/min state height' 'track variance' 'stasi states' 'stasi MDL' 'track number' 'track noise' 'track Quality' 'deviation from 2 state fit'}


%% state checking functions
function [states, exmid] = State_numbers(y, ideal_2, ideal_3, err_multi)

%inputs:
%y = height trace
%ideal_2 = 2 state ideal trace
%ideal_3 = 3 state ideal trace
%err_multi = factor to increase error thresholds by
%outputs:
%[states(j,1), states(j,2)] = State_numbers(track_hwd{j}, s2{j}.ideal, s3{j}.ideal,err_mult);

yup = 1*(ideal_2 > 0.99*max(ideal_2));
ydown = 1*(ideal_2 < min(ideal_2)+0.01);

up_values =  y(yup>0.5);
down_values =  y(ydown>0.5);

up_sd_2 = err_multi*std(up_values,0,2);

down_sd_2 = err_multi*std(down_values,0,2);
up_state_2 = max(ideal_2);
down_state_2 = min(ideal_2);

yup_3 = 1*(ideal_3 > 0.99*max(ideal_3));
ymid_3 = 1*(ideal_3>(min(ideal_3)+0.01) & ideal_3<0.95*max(ideal_3)); 
ydown_3 = 1*(ideal_3 < min(ideal_3)+0.01);

up_values_3 =  y(yup_3>0.5);
mid_values_3 = y(ymid_3>0.5);
down_values_3 =  y(ydown_3>0.5);

up_state_3 = max(ideal_3);
mid_state_3 = mean(ideal_3(ideal_3>(min(ideal_3)+0.01) & ideal_3<0.95*max(ideal_3))); 
down_state_3 = min(ideal_3);

up_sd_3 = err_multi*std(up_values_3,0,2);
mid_sd_3 = err_multi*std(mid_values_3,0,2); 
down_sd_3 = err_multi*std(down_values_3,0,2);

%if down_sd_3 < 0.15
%    down_sd_3 = 0.15;
%end
%if mid_sd_3 < 0.15
%  mid_sd_3 = 0.15;
%end
%if up_sd_3 < 0.15
%  up_sd_3 = 0.15;
%end

    %mostly up
    if numel(up_values_3)> numel(down_values_3) && numel(up_values_3)> numel(mid_values_3)
        
         if mid_state_3 <(up_state_3 - up_sd_3)  && mid_state_3 >(down_state_3 + down_sd_3)
        states = 3;
        elseif down_state_2 <(up_state_2 - up_sd_2) 
         states = 2;
        elseif down_state_2 > (up_state_2 - up_sd_2) 
         states = 1;
        else
          states = 1;
         end
    end
    
   %mostly down
   if numel(down_values_3)> numel(up_values_3) && numel(down_values_3)> numel(mid_values_3)
        
        if mid_state_3 >(down_state_3 + down_sd_3) && mid_state_3 <(up_state_3 + up_sd_3)
        states = 3;
        elseif up_state_2 >(down_state_2 + down_sd_2)
         states = 2;
        elseif up_state_2 <(down_state_2 + down_sd_2)
         states = 1;
        else
          states = 1;
        end
   end
   
   %mostly mid
   if numel(mid_values_3)> numel(up_values_3) && numel(mid_values_3)> numel(down_values_3)
        
        if down_state_3 <(mid_state_3 - mid_sd_3) && up_state_3 >(mid_state_3 + mid_sd_3)
        states = 3;
        elseif up_state_3 >(mid_state_3 + mid_sd_3)
         states = 2;
        elseif down_state_3 <(mid_state_3 - mid_sd_3)
         states = 2;
        else
         states = 1;
        end
      exmid = 0;
     else
      exmid = 1; 
   end
   exmid = states;%exmid*states-(exmid-1);
end
function [up_h, up_wd, mid_h, mid_wd, down_h, down_wd, ideal, ideal_i] = ThreeState_measure(x,y, ideal)

%x = track_time{3}; y = track_hwd{3};  ideal = output(3).eff_fit(3,:);
if any(y(:))
ideal_i = ideal;

yup = 1*(ideal > max(ideal)-0.01);
ymid = 1*(ideal>(min(ideal)+0.01) & ideal<(max(ideal)-0.01)); 
ydown = 1*(ideal < min(ideal)+0.01);

up_values =  y(yup>0.5);
mid_values = y(ymid>0.5);
down_values =  y(ydown>0.5);

up_state = max(ideal);
mid_state = mean(ideal(ideal>(min(ideal)+0.01) & ideal<0.95*max(ideal))); 
down_state = min(ideal);

up_sd = std(up_values,0,2);
mid_sd = std(mid_values,0,2); 
down_sd = std(down_values,0,2);

    [i_up_h,up_loc,up_wd] = findpeaks(yup,x,'MinPeakHeight',0.99,'MinPeakDistance',1,'MinPeakWidth',0.5);
    [i_mid_h,mid_loc,mid_wd] = findpeaks(ymid ,x,'MinPeakHeight',0.99,'MinPeakDistance',1,'MinPeakWidth',0.5);
    [i_down_h,down_loc,down_wd] = findpeaks(ydown ,x,'MinPeakHeight',0.99,'MinPeakDistance',1,'MinPeakWidth',0.5);

    for i=1:numel(mid_wd)
    pos = (x>(mid_loc(i)) & x<(mid_loc(i)+mid_wd(i)));
    mid_h(i,:) =  mean(y(pos));
    pos = (x>(mid_loc(i)-0.1) & x<(mid_loc(i)+mid_wd(i)));
        if mid_h(i,:)< mid_state
             if  mid_h(i,:)> (down_state+ down_sd) && (mid_state- mid_h(i,:))< ( mid_h(i,:)- down_state) && mid_wd(i)>1
             ideal(pos) =  mid_state;
             else
             ideal(pos) = down_state;
             end
        else
            if  mid_h(i,:)< (up_state- up_sd) && (mid_h(i,:)-mid_state)<(up_state- mid_h(i,:)) &&mid_wd(i)>1
            ideal(pos) =  mid_state;
            else
            ideal(pos) = up_state;
            end   
        end
    end 

    for i=1:numel(up_wd)
         pos = (x>(up_loc(i)-0.1) & x<(up_loc(i)+up_wd(i)));
         if up_wd(i) > 1
         ideal(pos) = up_state;    
         else
         ideal(pos) = down_state;
         end
    end
    
    for i=1:numel(down_wd)
         pos = (x>(down_loc(i)-0.1) & x<(down_loc(i)+down_wd(i)));
         if down_wd(i) > 1
         ideal(pos) = down_state;    
         else
         ideal(pos) = up_state;
         end
    end
    
mid_h = []; mid_wd = []; up_wd = []; down_wd = []; 
        
yup = 1*(ideal > max(ideal)-0.01);
ymid = 1*(ideal>(min(ideal)+0.01) & ideal<(max(ideal)-0.01)); 
ydown = 1*(ideal < min(ideal)+0.01);

    
%peak measuring
    [i_up_h,up_loc,up_wd] = findpeaks(yup,x,'MinPeakHeight',0.99,'MinPeakDistance',1,'MinPeakWidth',0.5);
    [i_mid_h,mid_loc,mid_wd] = findpeaks(ymid ,x,'MinPeakHeight',0.99,'MinPeakDistance',1,'MinPeakWidth',0.5);
    [i_down_h,down_loc,down_wd] = findpeaks(ydown ,x,'MinPeakHeight',0.99,'MinPeakDistance',1,'MinPeakWidth',0.5);

        %averaging heights in peak windows   
           if numel(up_wd) > 0
             for i=1:numel(up_wd) 
             pos = (x>(up_loc(i)-1) & x<(up_loc(i)+up_wd(i)-0.5));
             up_h(i,:) =  mean(y(pos)- min(ideal));     
             end  %max peaks
           else
            up_h = 0;
            end
        if numel(mid_wd)  > 0 
        for i=1:numel(mid_wd)  
        pos = (x>(mid_loc(i)-1) & x<(mid_loc(i)+mid_wd(i)-0.5));
        mid_h(i,:) =  mean(y(pos)- min(ideal));
        end  %mid peaks
       else
            mid_h = 0;
        end
            
       if numel(down_wd) > 0
        for i=1:numel(down_wd)
        pos = (x>(down_loc(i)-1) & x<(down_loc(i)+down_wd(i)-0.5));
        down_h(i,:) =  mean(y(pos)- min(ideal));
        end  %down peaks
       else 
           down_h = 0;
       end
end
end
function [ideal] = ThreeState_ToTwo(x,y, ideal)


if any(y(:))

yup = 1*(ideal > max(ideal)-0.01);
ymid = 1*(ideal>(min(ideal)+0.01) & ideal<(max(ideal)-0.01)); 
ydown = 1*(ideal < min(ideal)+0.01);


up_values =  y(yup>0.5);
mid_values = y(ymid>0.5);
down_values =  y(ydown>0.5);

up_state = max(ideal);
mid_state = mean(ideal(ideal>(min(ideal)+0.01) & ideal<0.95*max(ideal))); 
down_state = min(ideal);

up_sd = std(up_values,0,2);
mid_sd = std(mid_values,0,2); 
down_sd = std(down_values,0,2);
  if ymid(1) == 1
       x = vertcat(x(1)-1, x);
       ymid =[0 ymid];
        y =[0 y];
        ideal =[0 ideal];
        add_s=1;
         else 
        add_s=0;
  end
    if ymid(end) == 1
       x = vertcat(x, (x(end)+1));
       ymid =[ymid 0];
       y =[y 0];
       ideal =[ideal 0];
       add_e = 1;
    else 
        add_e=0;
  end
  
 [i_mid_h,mid_loc,mid_wd] = findpeaks(ymid,x,'MinPeakHeight',0.99,'MinPeakDistance',1,'MinPeakWidth',1.1);

    for i=1:numel(mid_wd)
    pos = (x>(mid_loc(i)) & x<(mid_loc(i)+mid_wd(i)));
    mid_h(i,:) =  mean(y(pos));
    pos = (x>(mid_loc(i)-0.1) & x<(mid_loc(i)+mid_wd(i)));
    
        if mid_h(i,:) - down_state < up_state - mid_h(i,:)
             ideal(pos) = down_state;
        else
            ideal(pos) = up_state;
        end
    end
if add_s ==1
    ideal(1) =[];
end
if add_e ==1
    ideal(end) =[];
end

end
end
function [up_h, up_wd, down_h, down_wd, ideal, ideal_i] = TwoState_measure(x,y, ideal)

%inputs:
%x = x of raw data
%y = heigth of raw data
%ideal = idealized 2 state trace of y

%script checks for errors in 2 state model for each state checking if it
%lies with the error of or closer to the other state based on the average 
%height of the raw data in that state 

%outputs:
%widths of states, average heights of states
%initail idealized trace = ideal_i
%treated idealized trace = ideal
if any(y(:))
ideal_i = ideal;
yup = 1*(ideal > (max(ideal)-0.01));
ydown = 1*(ideal < min(ideal)+0.01);

up_values =  y(yup>0.5);
down_values =  y(ydown>0.5);
up_sd = std(up_values,0,2);
down_sd = std(down_values,0,2);

%peak measuring
    [i_up_h,up_loc,up_wd] = findpeaks(yup,x,'MinPeakHeight',0.99,'MinPeakDistance',1,'MinPeakWidth',0.5);
    [i_down_h,down_loc,down_wd] = findpeaks(ydown ,x,'MinPeakHeight',0.99,'MinPeakDistance',1,'MinPeakWidth',0.5);

        for i=1:numel(up_wd)  
        pos = (x>(up_loc(i)-1) & x<(up_loc(i)+up_wd(i)-0.5));
        up_h(i,:) =  mean(y(pos));
        pos = (x>(up_loc(i)-0.1) & x<(up_loc(i)+up_wd(i)));
        
                   if  up_h(i,:)>(min(ideal)+ down_sd) && (up_h(i,:)-min(ideal))> (max(ideal) - up_h(i,:)) && up_wd(i)>1
                   ideal(pos) =  max(ideal);
                   else
                   ideal(pos) = min(ideal);
                   end
        end
        
        for i=1:numel(down_wd)
        pos = (x>(down_loc(i)-1) & x<(down_loc(i)+down_wd(i)-0.5));
        down_h(i,:) =  mean(y(pos));
        pos = (x>(down_loc(i)-0.1) & x<(down_loc(i)+down_wd(i)));
        
                   if  down_h(i,:)<(max(ideal)- up_sd) && (down_h(i,:)-min(ideal))< (max(ideal) - down_h(i,:))&& down_wd(i)>1
                   ideal(pos) =  min(ideal);
                   else
                   ideal(pos) = max(ideal);
                   end
        end
        
        up_h = []; up_wd = []; down_h=[]; down_wd=[];
        
yup = 1*(ideal > (max(ideal)-0.01));
ydown = 1*(ideal < min(ideal)+0.01);

        
        [i_up_h,up_loc,up_wd] = findpeaks(yup,x,'MinPeakHeight',0.99,'MinPeakDistance',1,'MinPeakWidth',0.5);
        [i_down_h,down_loc,down_wd] = findpeaks(ydown ,x,'MinPeakHeight',0.99,'MinPeakDistance',1,'MinPeakWidth',0.5);
        
        %averaging heights in peak windows   
        for i=1:numel(up_wd)  
        pos = (x>(up_loc(i)-1) & x<(up_loc(i)+up_wd(i)-0.5));
        up_h(i,:) =  mean(y(pos)- min(ideal));
        end  %max peak
        
        for i=1:numel(down_wd)
        pos = (x>(down_loc(i)-1) & x<(down_loc(i)+down_wd(i)-0.5));
        down_h(i,:) =  mean(y(pos)- min(ideal));
        end  %down peaks
end
end

function [ha, pos] = tight_subplot(Nh, Nw, gap, marg_h, marg_w)

% tight_subplot creates "subplot" axes with adjustable gaps and margins
%
% [ha, pos] = tight_subplot(Nh, Nw, gap, marg_h, marg_w)
%
%   in:  Nh      number of axes in hight (vertical direction)
%        Nw      number of axes in width (horizontaldirection)
%        gap     gaps between the axes in normalized units (0...1)
%                   or [gap_h gap_w] for different gaps in height and width 
%        marg_h  margins in height in normalized units (0...1)
%                   or [lower upper] for different lower and upper margins 
%        marg_w  margins in width in normalized units (0...1)
%                   or [left right] for different left and right margins 
%
%  out:  ha     array of handles of the axes objects
%                   starting from upper left corner, going row-wise as in
%                   subplot
%        pos    positions of the axes objects
%
%  Example: ha = tight_subplot(3,2,[.01 .03],[.1 .01],[.01 .01])
%           for ii = 1:6; axes(ha(ii)); plot(randn(10,ii)); end
%           set(ha(1:4),'XTickLabel',''); set(ha,'YTickLabel','')

% Pekka Kumpulainen 21.5.2012   @tut.fi
% Tampere University of Technology / Automation Science and Engineering


if nargin<3; gap = .02; end
if nargin<4 || isempty(marg_h); marg_h = .05; end
if nargin<5; marg_w = .05; end

if numel(gap)==1; 
    gap = [gap gap];
end
if numel(marg_w)==1; 
    marg_w = [marg_w marg_w];
end
if numel(marg_h)==1; 
    marg_h = [marg_h marg_h];
end

axh = (1-sum(marg_h)-(Nh-1)*gap(1))/Nh; 
axw = (1-sum(marg_w)-(Nw-1)*gap(2))/Nw;

py = 1-marg_h(2)-axh; 

% ha = zeros(Nh*Nw,1);
ii = 0;
for ih = 1:Nh
    px = marg_w(1);
    
    for ix = 1:Nw
        ii = ii+1;
        ha(ii) = axes('Units','normalized', ...
            'Position',[px py axw axh], ...
            'XTickLabel','', ...
            'YTickLabel','');
        px = px+axw+gap(2);
    end
    py = py-axh-gap(1);
end
if nargout > 1
    pos = get(ha,'Position');
end
ha = ha(:);
end


%% StaSI functions
%GNU General Public License v3.0
% see DOI: 10.1021/jz501435p
%https://github.com/LandesLab/STaSI

function output = StaSI_inputeff(eff, sd_control)
%% The main function of the generalized change point algorithm
% Input: 
%       single 1D trace (eff) or multiple traces selected in dialog
% Output:
%       G: structure, recording all the optimum clusterings under different
%   number of states
%       MDL: the minimum discription length, used to determine the optimum
%   number of states
%       states: the fitting based on the optimum number of states
%       eff: group the traces of different data together, if use input eff
%   trace, the output will be identical to the input
%       eff_fit: the fitting of all the feasible number of states, up to 30
%       breaks: recording the separations among different traces
%       output: recording several important parameters for potential usage
%       records: structure, recording the analysis of each loaded trace,
%   and also recording the location of each trace in the output eff
%       excluded: structure, recording all the traces not being used

%% step 1: loading the traces and change-points detection
sd = w1_noise(diff(eff))/1.4;% estimate the noise level
points = change_point_detection(eff);% change points detection
T2 = numel(eff);
breaks = T2;
groups = [1, points+1; points, T2];
try
    sd = sd_control;
catch ME
    sd = max(sd);% use the maximum noise level among these traces as the global noise level
end

%% step 2 and 3: clustering the segments and calculate MDL
[G, Ij, Tj] = clustering_GCP(eff, groups);
G = G(end:-1:1);% flip the G
n_mdl = min(30, numel(G));% calculate up to 30 states
MDL = zeros(1,n_mdl);
eff_fit = zeros(n_mdl, numel(eff));
for i = 1:n_mdl;
    [MDL(i), eff_fit(i,:)] = MDL_piecewise(Ij, Tj, G(i), eff, groups, sd, breaks);
end
%[~, q] = min(MDL);%now the BIC is actually MDL

output.MDL = MDL;
output.eff = eff;
output.eff_fit = eff_fit;
output.breaks = breaks;


end
%subfunctions
function sd = w1_noise(w1)
%% Estimate the standard deviation of noise from wavelet at scale = 1
% assum noise is Gaussian distribution, std of noise can be best estimated
% from w1 by counting.
% x ~ N(0, sig), then median of x is close to 0, and 34.1% away from median
% is about sig
% for abs(x), 68.2% away from the min is about the sig
y = abs(w1);
y = sort(y);
sd = y(round(0.682*numel(w1)));
% even for poisson distribution, as long as lambda >= 1, the median is
% almost = lambda. So, still can use counting method.
end
function points = change_point_detection(eff)
%% The main function to detect all the change points in a trace
sd = w1_noise(diff(eff))/1.4;% estimate the noise level
points = recursion1(eff,sd,[], 0);% recursively detect all the change points
points = sort(points);
end
function points = recursion1(eff, sd, points, counter)
%{
tau998_table = [1:30, 40, 50, 60, 80, 100, 1000;...
    318.31, 22.327, 10.215, 7.173, 5.893, 5.208, 4.785, 4.501, 4.297, 4.144,...
    4.025, 3.93, 3.852, 3.787, 3.733, 3.686, 3.646, 3.61, 3.579, 3.552,...
    3.527, 3.505, 3.485, 3.467, 3.45, 3.435, 3.421, 3.408, 3.396, 3.385,...
    3.307, 3.261, 3.232, 3.195, 3.174, 3.098];% the t-distribution
%}
N = numel(eff);
tau998 = 3.174;
if N < 2% only one point left in the segment, stop searching for the change point
    return
else
    llr = change_point_wavelet(eff,sd);
    [Z, k] = max(abs(llr));
    if Z > tau998
        counter = 0;
        points(end+1) = k;
        points1 = recursion1(eff(1:k), sd, [], counter);
        points2 = recursion1(eff(k+1:end), sd, [], counter);
        points = [points, points1, points2+k];
    elseif counter < 3% the parameter to dig in and find more short-lived transitions
        counter = counter +1;
        k = floor(numel(eff)/2);
        points1 = recursion1(eff(1:k), sd, [], counter);
        points2 = recursion1(eff(k+1:end), sd, [], counter);
        points = [points, points1, points2+k];
    else
        counter = 0;
        return
    end
end
end
function llr = change_point_wavelet(eff, sd)
N = numel(eff);
llr = zeros(size(eff));
for i = 1 : N-1
    I1 = mean(eff(1:i));
    I2 = mean(eff(i+1:end));
    llr(i) = (I2 - I1)/sd/sqrt(1/i+1/(N-i));
end
end
function [MDL, eff_fit] = MDL_piecewise(Ij, Tj, G, eff, groups, sd, breaks)
[Pmj, Im] = pre_calculation(Ij, Tj, G);
[nG, Ncp] = size(Pmj);% nG: the number of groups; Ncp: the number of change points detected
N = numel(eff);
V = max(eff)-min(eff);% the space of state values
eff_fit = zeros(size(eff));% the fitting
%nk = zeros(1,nG);
for k = 1:Ncp
    %nk(Pmj(:,k)==1) = nk(Pmj(:,k)==1)+groups(2,k)-groups(1,k)+1;
    eff_fit(groups(1,k):groups(2,k)) = Im*Pmj(:,k);
end
F = N*log(sd*sqrt(2*pi)) + 1/2/sd*sum(abs(eff-eff_fit));% the goodness of the fit
[lnDET, nb] = the_matrix(eff_fit, breaks, sd);
G = nG/2*log(1/2/pi)+nG*log(V/sd)+nb/2*log(N)+0.5*lnDET;% the cost of the model
MDL = F+G;
end
function [lnDET, nb] = the_matrix(states, breaks, sd)
[Z, ~, IZ] = unique(states);
nz = numel(Z);
indx = find(diff(IZ)~=0);%indx records the positions of transitions, and the period of each stay
for i = 1:numel(breaks)
    indx(indx == breaks(i)) = [];%these breaks will not be considered as a change point
end
% lucky, the matrix is diagonal
nb = numel(indx);
A = (zeros(1,nz));% the state value part
D = (zeros(1,nb));% the change-points part
for i = 1:nz
    A(i) = numel(IZ(IZ==i));% length of each state
end
for i = 1:nb
    temp = states(indx(i))-states(indx(i)+1);% the transition value
    D(i) = temp^2;
end
lnDET = sum(log([A,D/sd^2]));% let D divide sd here to avoid super large/small numbers
end
function [Pmj, Im] = pre_calculation(Ij, Tj, G)
N = numel(Ij);
n = numel(G.g);
Pmj = zeros(n, N);
for i = 1:n
    for j = 1:numel(G.g(i).gg)
        Pmj(i,G.g(i).gg(j))=1;
    end
end
Tm = zeros(1,n);
nm = zeros(1,n);
Im = zeros(1,n);
for i = 1:n
    Tm(i) = sum(Pmj(i,:).*Tj);
    nm(i) = sum(Pmj(i,:).*(Ij.*Tj));
    if Tm(i) ~= 0
        Im(i) = nm(i)/Tm(i);
    else
        Im(i) = 0;
    end
end
end
function [G, Ij, mj, split_tree] = clustering_GCP(eff, groups)
%% clusetering all the segments one by one up to only one
len = length(groups(1,:));% number of segments
mj = groups(2,:)-groups(1,:)+1;% the number of data points in each segment
Ij = zeros(1,len);% the averaged intensity of each segment
G = struct([]);% the structure to store all the clustering history
for i = 1:len
    Ij(i) = mean(eff(groups(1,i):groups(2,i)));
    G(1).g(i).gg = i;
end
[G, split_tree] = sub_clustering(Ij, mj, G);% clustering
end
function [G, split_tree] = sub_clustering(Ij, mj, G)
n = numel(G(end).g);
if n <= 1
    split_tree = [];
    return
else
    M = ones(n, n)*(-inf);% the merit matrix
    for i = 1:n-1
        for j = i+1:n
            % calculating the merit for all the possible pairs
            M(i, j) = llr_merit(Ij(i), Ij(j), mj(i), mj(j));
        end
    end
    temp = [];
    while numel(G(end).g) > 1;
        numel(G(end).g)
        [a, b] = size(M);
        [~, q] = max(M(:));
        j = ceil(q/a);% find the column
        i = q-(j-1)*a;% find the row
        %%
        if numel(G(end).g) <= 30%% to record the split tree
            n = numel(G(end).g);
            if isempty(temp)
                temp = ones(1, n)*(n+0.5);
                split_tree = zeros(4, 3*n-1);
                split_tree(1:2, 2*n:end) = n+0.5;
                split_tree(3, 2*n:end) = Ij;
                split_tree(4, 2*n:end) = mj;
            end
            split_tree(1, 2*n-2:2*n-1) = n;
            split_tree(2, 2*n-2:2*n-1) = temp([i,j]);
            split_tree(3, 2*n-2:2*n-1) = Ij([i,j]);
            split_tree(4, 2*n-2:2*n-1) = mj([i,j]);
            temp(i) = n;
            temp(j) = [];
        end
        G(end+1).g = G(end).g;% initialize G(n+1)
        G(end).g(i).gg = [G(end).g(i).gg, G(end).g(j).gg];% group segments i and j together
        G(end).g(j) = [];% free another space
        Ij(i) = (Ij(i)*mj(i)+Ij(j)*mj(j))/(mj(i)+mj(j));% group Ij
        Ij(j) = [];
        mj(i) = mj(i)+mj(j);% group mj
        mj(j) = [];
        %% update the M
        M(:, i) = -inf; M(i, :) = -inf;
        M(:, j) = []; M(j, :) = [];
        for k = 1:numel(G(end).g)% pairing with the new cluster i
            if k < i
                M(k, i) = llr_merit(Ij(i), Ij(k), mj(i), mj(k));
            elseif k > i
                M(i, k) = llr_merit(Ij(i), Ij(k), mj(i), mj(k));
            end% if
        end% for k
    end% while
    split_tree(1:4,1) = [1; 2; Ij; mj];
end% if n<=1
end% function
function llr = llr_merit(I1, I2, m1, m2)
I = (I1*m1 + I2*m2)/(m1+m2);
llr = (m1+m2)*I^2-m1*I1^2-m2*I2^2;
end

