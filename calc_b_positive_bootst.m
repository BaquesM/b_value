function [b_value,b_value_std]=calc_b_positive_bootst(cat,nb_bootst,mag_cor)
% INPUT
% cat: must contain magnitude
% nb_bootst: number of iteration bootstrapping
% mag_cor: Linked to the difference of magnitude completness obtained by MAXC; Generally fixed at 0.2
% OUTPUT
% b_value: b-value calculated from b-positive method (van der Elst (2021))
% b_value_std: difference between the b-value calculated and the mean of the b-value calculated with bootstrap

% COMMENTS
% The catalogue has to be sort by time

% Created by: Marion BAQUES
% Version: 17/08/2022

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%																                                                                 %
%				                    CALCULATION OF B-POSITIVE VALUE FOR ALL CATALOGUES						                     %
%																                                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculation of the difference between two magnitudes
cpt=0;
for i=2:length(cat)
    cpt=cpt+1;
    Diff_mag(cpt)=cat(i)-cat(i-1);
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

% Positive magnitude difference
Diff_mag_pos=Diff_mag(Ind_m);

if length(Ind_m) >= 70 
    if nb_bootst < 2
        % Calculation of the beta positive estimator
        Beta=(mean(Diff_mag_pos)-Min_pos_mag_diff)^-1;
        % Calculation of b-value
        b_value=Beta/log(10);
        % Standard error for the b-value
        b_value_std=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%															                                                                	 %
%				                CALCULATION OF UNCERTAINTY B-POSITIVE VALUE FOR ALL CATALOGUES					                 %
%															                                                                	 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    else
        for nb=1:nb_bootst
	        clear min_diff_pos max_diff_pos Diff_perturb
            	disp(['nb bootstrap: ',num2str(nb),'/',num2str(nb_bootst)])
	        min_diff_pos=min(Diff_mag_pos);
	        max_diff_pos=max(Diff_mag);
	        Diff_perturb=Diff_mag_pos(randi(length(Diff_mag_pos),length(Diff_mag_pos),1));
		% Calculation of the beta positive estimator
		Beta(nb)=(mean(Diff_perturb)-Min_pos_mag_diff)^-1;
		% Calculation of b-value
		b_pos(nb)=Beta(nb)/log(10);
	    end
        % Standard error for the b-value
        b_value=mean(b_pos);
        b_value_std=std(b_pos)*2;
    end    
    % % Figure b_value
    %edges=[b_value-0.4:0.01:b_value+0.4];
    %histogram(b_pos,edges)
    %xlabel('Perturbed b-value')
    %ylabel('Nb of b-value')
    %box on
    %set(gca,'FontSize',22)
    %saveas(1,'Figure_b_value_distribution','pdf')
else
    disp(['Insuffisant number of events: ', num2str(length(Ind_m))])
    b_value=0;
    b_value_std=0;
end


