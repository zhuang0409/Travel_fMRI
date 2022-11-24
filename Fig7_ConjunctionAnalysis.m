%% Conjunction maps
% to identify regions in which actions are jointly represented (a) at the 
% subordinate and basic level, and (b) at all three taxonomic levels, we 
% computed a conjunction analysis on the basis of the statistical maps of 
% RSA results with regressing out the visual model (the corresponding 
% results shown in Figure 6A). 

%% FOLDER
projectdir='G:\Travel\data\ExemData\RSA\RSA_glmRN50\Group\';
mask_fn = 'G:\Travel\data\ExemData\MNI152_T1_2mm_brain_mask.nii.gz';

%% read data of each level 
% maps after TFCE
% super
super_dsfn=fullfile(projectdir,sprintf('TFCE_super_RN50_statmap.nii.gz'));
ds_super = cosmo_fmri_dataset(super_dsfn,'mask',mask_fn);
% basic
basic_dsfn=fullfile(projectdir,sprintf('TFCE_basic_RN50_statmap.nii.gz'));
ds_basic = cosmo_fmri_dataset(basic_dsfn,'mask',mask_fn);
% sub
sub_dsfn=fullfile(projectdir,sprintf('TFCE_sub_RN50_statmap.nii.gz'));
ds_sub = cosmo_fmri_dataset(sub_dsfn,'mask',mask_fn);
%% conjunction
% Conjunctions were computed by outputting the minimum t value for each vertex of the input maps
% conjunct of sub and basic level: to get the min t between basic and sub
ds_conj_SubBasic=min(ds_basic.samples,ds_sub.samples);
% conjunct of three taxonomic level: to get the min t across super, basic, sub.
ds_conj=min(ds_basic.samples,ds_super.samples);
ds_conj_all=min(ds_conj,ds_sub.samples);
%% save maps
super_map.samples=ds_conj_all;
super_map.a=ds_super.a;
super_map.sa=ds_super.sa;
super_map.fa=ds_super.fa;

basic_map.samples=ds_conj_SubBasic;
basic_map.a=ds_basic.a;
basic_map.sa=ds_basic.sa;
basic_map.fa=ds_basic.fa;

cosmo_map2fmri(super_map, fullfile(projectdir, sprintf('GROUP_SuperBasicSub.nii.gz')));
cosmo_map2fmri(basic_map, fullfile(projectdir, sprintf('GROUP_BasicSub.nii.gz')));