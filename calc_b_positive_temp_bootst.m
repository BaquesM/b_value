function [b_value_temp,b_value_temp_std,Time_beg,Time_end]=calc_b_positive_temp_bootst(cat,mag_cor,nb_bootst,bin,ini,dt)
% INPUT
% cat: must contain magnitude
% mag_cor: Linked to the difference of magnitude completness obtained by MAXC; Generally fixed at 0.2
% nb_bootst: number of iteration bootstrapping
% Input format of cat: YYYY MO DY HR MN SC LON LAT DEP MAG
% Input format of Unc: ID YYYY MO DY HR MN SC LON LAT DEP MAG
% bin: bin of positive magnitude difference
% ini: number of positive magnitude difference in commom between bin
% dt: number of events not in commom between bin
% OUTPUT
% b_value: b-value calculated from b-positive method (van der Elst (2021))
% b_value_std: difference between the b-value calculated and the mean of the b-value calculated with bootstrap
% Time_beg: Starting time of the window
% Time_end: Ending time of the window

% COMMENTS
% The catalogue has to be sort by time

% Created by: Marion BAQUES
% Version: 22/08/2022

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%																                                                                 %
%				                    CALCULATION OF TEMPORAL B-POSITIVE VALUES FOR CATALOGUE						                 %
%																                                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cpt=0;
% Calculation of the magnitude difference
for i=2:length(cat(:,1))
    cpt=cpt+1;
    Diff_mag(cpt)=cat(i,10)-cat(i-1,10);
    % Time of the first and second events  (af JC)
    Time_border(cpt,:)=[datenum(cat(i-1,1),cat(i-1,2),cat(i-1,3),cat(i-1,4),cat(i-1,5),cat(i-1,6)) datenum(cat(i,1),cat(i,2),cat(i,3),cat(i,4),cat(i,5),cat(i,6))];
end
% If no uncertainty file available, use the general correction
if mag_cor==0
    mag_cor=0.2;
else
    mag_cor=mag_cor;
end

% Find the minimum positive magnitude difference
Ind=find(Diff_mag>0);
Min_pos_mag_diff=min(Diff_mag(Ind))+mag_cor;
% Find the difference higher or equal to the minimum positive magnitude difference
Ind_m=find(Diff_mag>=Min_pos_mag_diff);

% Time since JC
Time_cat=datenum(cat(:,1),cat(:,2),cat(:,3),cat(:,4),cat(:,5),cat(:,6));
% Positive magnitude difference
Diff_mag_pos=Diff_mag(Ind_m);

%%%%%%%%%%%%%%%%%%%%%%%%%%% Calculation of uncertainty %%%%%%%%%%%%%%%%%%%%%%%%%%% 
% For each bin of events
clear Same_evt Cat unc_mag
if nb_bootst > 1
    % For each bin of 50 positive Difference magnitude greater than the minimum positive magnitude, calculate the b-value 
    cpt3=0;
    for i=ini:dt:length(Ind_m)     
        if (i-(ini-1))+(bin-1) <= length(Ind_m)
	        cpt3=cpt3+1;
            %disp(['iter: ',num2str(cpt3)])
            for nb=1:nb_bootst
	            clear min_diff_pos max_diff_pos Diff_perturb
		        di=Diff_mag_pos(i-(ini-1):(i-(ini-1))+(bin-1));
	            min_diff_pos=min(di);
	            max_diff_pos=max(di);
                Diff_perturb=di(randi(length(di),length(di),1));
                % Calculation of the beta positive estimator
                Beta(nb)=(mean(Diff_perturb)-Min_pos_mag_diff)^-1;
                % Calculation of b-value
                b_pos(nb)=Beta(nb)/log(10);
	        end
	        b_value_temp(cpt3)=mean(b_pos);
	        b_value_temp_std(cpt3)=std(b_pos)*2;
	        Time_beg(cpt3)=Time_border(Ind_m(i-(ini-1)),1);
	        Time_end(cpt3)=Time_border(Ind_m((i-(ini-1))+(bin-1)),1);
        end
    end
else
    % Without uncertainty value of magnitude
    % For each bin of 30 positive Difference magnitude greater than the minimum positive magnitude, calculate the b-value 
    cpt3=0;
    for i=ini:dt:length(Ind_m)
        if (i-(ini-1))+(bin-1) <= length(Ind_m)
            cpt3=cpt3+1;
            clear Beta BB
            BB=Diff_mag(Ind_m(i-(ini-1):(i-(ini-1))+(bin-1)));
            % Calculation of the beta positive estimator
            Beta=(mean(BB)-Min_pos_mag_diff)^-1;
            % Calculation of b-value
            b_value_temp(cpt3)=Beta/log(10);
            b_value_temp_std(cpt3)=0;
            Time_beg(cpt3)=Time_border(Ind_m(i-(ini-1)),1);
            Time_end(cpt3)=Time_border(Ind_m((i-(ini-1))+(bin-1)),1);
        end
    end   
end
