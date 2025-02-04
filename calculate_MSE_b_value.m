function [MSE,MSE_percent,N_data_save,N_rebuild,v_rebuild]=calculate_MSE_b_value(b_value,Diff_pos_bin,nbin)
% INPUT:
% b_value: b-value of Gutenberg-Richter
% min_mag: minimum positive magnitude difference (from b-positive method, van der Elst (2021))
% max_mag: maximum positive magnitude difference (from b-positive method, van der Elst (2021))
% Diff_pos_bin: positive magnitude difference (for each calculated b-value)
% nbin: number of positive magnitude difference in a b-value
% OUTPUT:
% MSE: mean squared error
% MSE_percent: mean squared error in percentage
% N_data_save: positive magnitude difference used to build GR
% N_rebuild: positive magnitude difference reconstructed to build GR
% v_rebuild: positive magnitude difference range reconstructed to build GR

% Created by: Marion BAQUES
% Version: 04/02/2025

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%																                                             %
%				                    CALCULATION OF MSE FOR B-VALUE						                     %
%																                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialization of parameters
MSE=zeros(length(b_value(1,:)),1);
MSE_percent=zeros(length(b_value(1,:)),1);
for nb_b=1:length(b_value(1,:))
    % Cumulative Difference magnitude positive
    clear nZ cZ indx ind_diff N_data first_part second_part
    nZ = hist(Diff_pos_bin(nb_b, :), v);
    for i = 1:length(nZ)
        cZ(i) = sum(nZ(i:end));
    end
    % 90% max N and more than N=10
    ind_diff=find(cZ> 10 & cZ <= ((nbin*90)/100));
    N_data=cZ(ind_diff);
    N_data_save(nb_b,:)=cZ;
    a_value(nb_b,1)=mean(log10(N_data)+(b_value_temp(1,nb_b)*v(ind_diff)));
    % Calculation of N_rebuild for each Delta_mag
    N_rebuild(nb_b,1:size(N_data,2))=10.^(a_value(nb_b,1)-(b_value_temp(1,nb_b)*v(ind_diff)));
    v_rebuild(nb_b,1:size(N_data,2))=v(ind_diff);
    % Calculation MSE
    n_rebuild=N_rebuild(nb_b,:);
    n_rebuild(n_rebuild==0)=[];
    first_part=(N_data-n_rebuild).^2;
    first_part(first_part==0)=[];
    second_part=sum(first_part);
    MSE(nb_b,1)=(second_part/sum(N_data.^2));
    MSE_percent(nb_b,1)=(second_part/sum(N_data.^2))*100;
end
