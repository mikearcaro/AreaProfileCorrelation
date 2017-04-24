function [mean_within_profile, mean_crosshemi_roi_corr] = AreaProfileCorr
% AreaProfileCorr Area profile correlation analysis reported in
% Arcaro & Livingstone 2017, "A hierarchical, retinotopic proto-organization
% of the primate visual system at birth"
%
% [mean_within_profile, mean_crosshemi_roi_corr] = AreaProfileCorr plots and returns
% within and between hemisphere correlation matrices for all visual ROIS.
% The function also writes out 2 fMRI volumes (BRIK format):
% The files have the suffix _rh and _lh. Each contains whole-brain 
% areal profile correlations based on visual ROI seeds from the corresponding hemisphere.
% For _rh, the left hemisphere comprises the cross hemisphere correlations
% reported in the paper, and vice versa for _lh. 
% For each area, a correlation coefficient, t-stat, and areal profile 
% similarity coefficient are calculated. This latter coefficient corresponds to the data 
% reported in Figure 3. 
%
% The function loads in sample data, ROIs, and a brainmask. The sample
% data are a subset of the data reported in the paper, but will produce
% similar results. 
%
% The function requires the AFNI Matlab toolbox for the final step of
% writing out a BRIK volume. You can get the toolbox for free here:
% https://afni.nimh.nih.gov/
%
% Written by Mike Arcaro 2017


% These parameters can be edited - defaults set to paper analysis
out_filename=['baby4_areaprofile_example_subsample']; % filename for volume ouput
roi_iterations=100;         % number of iterations for subsampling ROI set to 1 for no subsampling
sim_corrfunc='Pearson';     % set correlation type for last step of analysis - Pearson / Kendall / Spearman
leave_self_out=0;           % set to 1 if you want to compute profile correlations w/o self-correlations for within areal correlation profile


% Do not change parameters below 
ROI_labels={'V1';'V2';'V3';'V4';'V4A';'MT';'MST';'FST';'V4t';'OTd';'PITv';'PITd';'OTS';'V3A/DP';'CIP';'LIP'}; % Area labels
hemis={'rh','lh'}; % define hemisphere order - needs 2 hemispheres

% Load data, brainmask, rois
disp('Loading data, brainmask, ROIs...')
load('brainmask')       % load brainmask - only ocmpute correlation within brain
load('ROIs')            % load rois   
load('data')            % load demo data

% Initialize variables
num_rois=size(ROI_labels,1);
data_length=size(data,2);

% Find minimum ROI size for each hemisphere and use to subsample
% the voxel selection of each ROI for correlation analysis
if roi_iterations>1
    roi_tc=zeros(roi_iterations,size(data,2),num_rois,size(hemis,1));
    % find minimum ROI size within hemisphere and use for voxel subsampling
    for h = 1:size(hemis,2)
        min_roi_size(h)=min(histc(ROIs{h},[1:num_rois])); % get voxel count for each ROI and find minimum
    end
else
    roi_tc=zeros(size(data,2),num_rois,size(hemis,1));
end

%% Compute crosshemisphere temporal correlations
% Equivalent to standard seed-based correlation anlaysis
% Seeds are individual visual areas correlated with whole brain
disp('Running Temporal Correlations...')
for h = 1 : size(hemis,2)
    % initialize temporary variables
    if roi_iterations>1
        crosshemi_corr_iter=zeros(size(data,1),roi_iterations,num_rois);
    else
        crosshemi_corr_temp=zeros(size(data,1),num_rois);
    end
    
    % loop through each ROI
    for r = 1:num_rois        
        if roi_iterations>1
            for curr_iter = 1:roi_iterations
                curr_roi=find(ROIs{h}==r); % select current ROI indices
                permnums=randperm(size(curr_roi,1)); % permute indices
                curr_roi_sample=curr_roi(permnums(1:min_roi_size(h))); % select first 1:minROIsize of permuted indices
                mask_index=ismember([1:size(data,1)],curr_roi_sample); % select data rows corresponding to ROI indices
                roi_tc(curr_iter,:,r,h)=mean(data(mask_index,:)); % calculate mean signal from those ROI indices
            end
            crosshemi_corr_iter(:,:,r)=corr(data',roi_tc(:,:,r,h)'); % correlate the mean ROI signal w/ whole brain data
        else
            roi_tc(r,:,h)=mean(data(ROIs{h}==r,:));
            crosshemi_corr_temp(:,r)=corr(data',roi_tc(r,:,h)');
        end    
    end
    
    % get average across iterations to get mean cross-hemi correlation
    % if iter > 1, average across subsampled correlations
    if roi_iterations>1
        crosshemi_corr(:,:,h)=squeeze(mean(crosshemi_corr_iter,2));
    else
        crosshemi_corr(:,:,h)=crosshemi_corr_temp;
    end
end

%% Compute the within-hemisphere areal correlation profile temporal correlations between all area pairs
disp('Computing Area Profile Similarity...')
% initialize temporary variables
if roi_iterations>1
    within_profile_iter=zeros(num_rois,num_rois,size(hemis,2),roi_iterations);
end

for h = 1:size(hemis,2)
    if roi_iterations>1
        for curr_iter = 1:roi_iterations
            within_profile_iter(:,:,h,curr_iter)=corr(squeeze(roi_tc(curr_iter,:,:,h)),squeeze(roi_tc(curr_iter,:,:,h)));
        end
        within_profile=mean(within_profile_iter,4);
    else
        within_profile(:,:,h)=corr(roi_tc(:,:,h),roi_tc(:,:,h));
    end
end
mean_within_profile=mean(within_profile,3);

figure 
imagesc(mean_within_profile)
title('Within Hemisphere Area Correlations')


%% Correlate the within hemisphere areal correlation profile with the voxel-wise cross-hemisphere areal correlation profile
for h = 1:size(hemis,2)
    if leave_self_out
        within_profile_temp=within_profile(:,:,h)'; within_profile_temp(logical(eye(size(within_profile_temp))))=NaN;
        crosshemi_sim(:,:,h)=corr(crosshemi_corr(:,:,h)',within_profile_temp,'type',num2str(sim_corrfunc),'rows','pairwise');
    else
        crosshemi_sim(:,:,h)=corr(crosshemi_corr(:,:,h)',within_profile(:,:,h),'type',num2str(sim_corrfunc));
    end
end

%% Create cross hemisphere correlation matrices
if roi_iterations>1
    for curr_iter = 1:roi_iterations
        crosshemi_roi_corr_iter(:,:,curr_iter)=corr(squeeze(roi_tc(curr_iter,:,:,1)),squeeze(roi_tc(curr_iter,:,:,2)));
        withinhemi_crosscorr_iter(:,:,curr_iter)=corr(within_profile_iter(:,:,1,curr_iter),within_profile_iter(:,:,2,curr_iter));
    end
    mean_crosshemi_roi_corr=mean(crosshemi_roi_corr_iter,3);
    mean_withinhemi_crosscorr=mean(withinhemi_crosscorr_iter,3);
else
    mean_crosshemi_roi_corr=corr(roi_tc(:,:,1),roi_tc(:,:,2));
    mean_withinhemi_crosscorr=corr(within_profile(:,:,1),within_profile(:,:,2));
end

figure
imagesc(mean_crosshemi_roi_corr)
title('Cross Hemisphere Temporal Correlations')



%% Combine seed correlations w/ stats for writing out to volume
corr_output=cell(2);
labels=cell(2);
for h = 1:size(hemis,2)
    for r = 1:num_rois
        crosshemi_corr_mean=crosshemi_corr(:,r,h);
        t_stats=(crosshemi_corr_mean)./sqrt( ( (1-crosshemi_corr_mean.^2)/(data_length-2)));
        corr_output{h}=[corr_output{h}, crosshemi_corr_mean, t_stats];
        labels{h}=[labels{h} num2str(ROI_labels{r}) '~ ' num2str(ROI_labels{r}) '_t~ '];
        
        crosshemi_sim_curr=crosshemi_sim(:,r,h);
        corr_output{h}=[corr_output{h}, crosshemi_sim_curr];
        labels{h}=[labels{h} num2str(ROI_labels{r}) '_sim~ '];
    end
end

%% Write out volume in BRIK format - requires AFNI matlab toolbox
disp('Writing data to BRIK...')
for h = 1:size(hemis,2)
    corr_output_full=zeros(size(brainmask.mask,1), size(corr_output{h},2));% analysis was on the subset of voxels that fell within brain
    corr_output_full(brainmask.mask==1,:)=corr_output{h}; % now we have to recreate the original volume based on the brainmask
    corr_output_vol=reshape(corr_output_full,[ brainmask.dims(1), brainmask.dims(2), brainmask.dims(3), size(corr_output_full,2)]);
    
    Info=brainmask.info; % use defaults from the brainmask
    Info.DATASET_DIMENSIONS(4)=size(corr_output_vol,4);
    Info.DATASET_RANK(2)=size(corr_output_vol,4);
    Info.BRICK_STATS=ones(1,size(corr_output_vol,4)*2);
    Info.BRICK_TYPES=3*ones(1,size(corr_output_vol,4));
    Info.BRICK_FLOAT_FACS=zeros(1,size(corr_output_vol,4));
    Info.TAXIS_NUMS(1)=Info.DATASET_RANK(2);
    Info.BRICK_LABS=labels{h};
    Opt.OverWrite='y';
    Opt.Prefix=[out_filename '_' num2str(hemis{h})];
    
    Opt.view ='orig';
    WriteBrik(corr_output_vol,Info,Opt);
end


end