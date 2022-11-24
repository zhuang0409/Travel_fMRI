% ROI-based RSA (partial correlation)
% set path
addpath('/Users/zhuang/EEGNet/results_ica/npy-matlab-master/npy-matlab');
%% FOLDER
projectDir = '/Users/zhuang/Documents/MRI/Projects/Travel/data/ExemData/RSA/ROIs_based/RSA_glm/';
mkdir (projectDir);
dataPath = '/Users/zhuang/Documents/MRI/Projects/Travel/data/ExemData/UnSmoothed/Exemdata/';
dsmdir='/Users/zhuang/Documents/MRI/Projects/Travel/data/ExemData/progs/';
%% Load ROI
roisdir='/Users/zhuang/Documents/MRI/Projects/Travel/data/ExemData/ROIs/10mm/';
rois=dir([roisdir, '*.nii.gz']);
smoothingVal=0;
% subjects
subjectIDs = {'SUB03_19980219SNFS','SUB04_19900101WALE','SUB05_19890101WANL','SUB06_19880720WAVI'...
              'SUB07_19960420WIST','SUB08_19980101THAE','SUB09_20200828NICA','SUB10_20200828LYXU'...
              'SUB11_19920409THZH','SUB12_19980908SABA','SUB13_19940216NARA','SUB14_19971002COCA'...
              'SUB15_19970428MIRU','SUB16_19891030CHZH','SUB17_19921010XIHA','SUB18_19921211ZUKA'...
              'SUB19_19970603JOBE','SUB20_19970125FIGI','SUB21_19940526MISC','SUB22_19891024ROPU'...
              'SUB23_19811010CHZW','SUB24_20200918ANIO','SUB25_20200923MICA'};

nsubjects=numel(subjectIDs);
levels=[{'super'} {'basic'} {'sub'}];

%% load the visual model
modeldir='/Users/zhuang/Documents/MRI/Projects/Travel/data/Resnet50/outputs/';
matrixnames=dir([modeldir,'*.npy']);
data= readNPY([modeldir,matrixnames(1).name]);
modeldsm=squareform(VGG_reshape_to_fMRI(squeeze(data(1,:,:))));

%% run partial correlation separately one by one
for i=1:length(rois)
    % load data
    mask_fn=fullfile([roisdir,rois(i).name]);
    display(rois(i).name);
    parital_r=zeros(nsubjects,length(levels));
    for j=1:nsubjects
        subjectID=subjectIDs{j};
        display(subjectID);
        data_fn=fullfile([dataPath, sprintf('Group_searchlight_SS%d_%s_aligned-subject-standard-2mm.nii.gz', smoothingVal, subjectID)]);
        ds_all=cosmo_fmri_dataset(data_fn,'mask',mask_fn);
        ds_all=cosmo_remove_useless_data(ds_all);
        img=repmat((1:72)',1,6);
        ds_all.sa.img=img(:);
        ds= cosmo_fx(ds_all,@(x)mean(x,1), 'img', 1);  
        
        ds_rdm=1-corr(ds.samples');
        ds_rdm_vec=squareform(ds_rdm);
        
        for k=1:3%super-1;basic-2;sub-3
            level=char(levels(k));
            load([dsmdir,level '.mat']);
            dsm_vec=squareform(dsm);
            X=[ds_rdm_vec' dsm_vec' modeldsm'];
            r=partialcorr(X,'type','Spearman');
            parital_r(j,k)=r(1,2);
        end
    end
        save ([projectDir, 'glmRSA_' rois(i).name(1:end-7) ], 'parital_r');
end
