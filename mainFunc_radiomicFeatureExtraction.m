function [features_extracted] = mainFunc_radiomicFeatureExtraction(workingDir, relative_path_to_ROI, outputDir, prefix_of_output_csv)
% Author: Chao Huang (huangchao312@gmail.com, huangchao09@zju.edu.cn)
% updated: May 17, 2018.
% Sourcefile: mainFile_feature_extraction_from_DICOM_n_DSO.m

% Inputs:
%   - workingDir: the working directory, where you store the images, DSO files and
%   want to store the results.
%   - relative_path_to_ROI: the path to ROI folder relative to the
%   workingDir. This folder only have all the images for the studied
%   patients, with images of each patient stored in a separate subfolder.

% workingDir =
% '/Users/messi/Documents/private_research/Stanford_BMI/TCIA_TCGA/';
cd (workingDir);

% relative_path_to_ROI = 'NSCLC Radiogenomics/';
% prefix_of_output_csv = 'NSCLC_radiogenomics_radiomicFeatures_';

% NOTE(1): Instead of using dicom_read_volume() written by the anonymous, the volume preparation code was adapted from
% <https://github.com/mvallieres/radiomics/>, except that double() was not
% applied. Some of the volumes are the same while some are not. Don't know
% why. Chao.

% NOTE(2): Within RaidomicsToolBox_Chao collection, I noticed getGaborBank() applied filter to getGlobalTextures(),
% N=3; getGLCMtextures(), N=11; getGLRLMtextures(), N=13;
% getGLSZMtextures(), N=13; getNGTDMtextures(), N=5, for each of 8 Frequencies and 8 Directions
% So final N = (8+8)*(3+9+13+13+5)=688. NOTICE: getGlobalTextures(),
% getGLCMtextures(), and getGLRLMtextures() are not the same functions to
% compute the corresponding feature groups without filtering. Not sure if I
% should make up for this. Chao.



%% initializers 
%filepath = '/Users/mzhou/Documents/temp_project/Tej/rois';
if exist('temp_file', 'dir')
    rmdir('temp_file');
end
mkdir('temp_file');

filepath = strcat(workingDir, relative_path_to_ROI);
All_Data_Path = dir(filepath);
% delete '.', '..', '.DS_Store', and '_DS_Store'
toRemoveIdx = [];
for i = 1:length(All_Data_Path)
    if strcmp(All_Data_Path(i).name,'.') || strcmp(All_Data_Path(i).name,'..') || strcmp(All_Data_Path(i).name,'.DS_Store') || strcmp(All_Data_Path(i).name,'_DS_Store') || strcmp(All_Data_Path(i).name,'._.DS_Store') % '._.DS_Store' required for Linux.
        toRemoveIdx = [toRemoveIdx,i];
    end;
end;
All_Data_Path(toRemoveIdx,:) = [];

formatOut = 'yyyymmdd'; % format date time to be like '20171225'

% end


%% collect dicom headers for QC, mask and volume for later use. The calculatation to calculate slice spacing is the only method available
% for now. Refer to: <https://github.com/mvallieres/radiomics/>
dicomInfo = struct;
% column 1, all slices' headers; column 2, computed sliceS(slice spacing).
% column 3, slice thickness.
for j =1:length(All_Data_Path) % modified by Chao
    fprintf('Reading headers, mask and volume from %d th case \n',j);
    %load data
    Directory = fullfile(filepath, All_Data_Path(j).name);
    elements=subdir(Directory); % list folders and subfolders. subdir()
    % Copyright (c) 2015 Kelly Kearney.
    nElements = length(elements); % includes all elements:e.g. .DS_Store files, mask file.
    volume = cell(1,1,nElements); % initialize the volume depth the same as the number of all elements. will trim later.
    dicomHeaders = []; % initialization, dicomHeaders refer to headers for
    % all DICOM files from a series
    SEG = []; % will store mask array from DSO file. The modality is "SEG".
    
    dicomInfo(j).patID = All_Data_Path(j).name; % patient ID
    
    % starts to read
    sliceNumber = 0;
    for elementNumber = 1:nElements
        elementName = elements(elementNumber).name;
        if ~strcmp(elementName,'.') && ~strcmp(elementName,'..') && ~contains(elementName,'.DS_Store') && ~contains(elementName,'_DS_Store') && ~contains(elementName,'._.DS_Store')
            % Good enough for Linux, add conditions for MAC and Windows.
            % strcmp() compares indentical strings. contains() means matching a
            % pattern.
            elementFullFile = elementName;
            if isdicom(elementFullFile)
                tmp = dicominfo(elementFullFile);
                if strcmp(tmp.Modality,'RTSTRUCT')
                    RTstruct = tmp;
                elseif strcmp(tmp.Modality,'REG')
                    REG = tmp;
                elseif strcmp(tmp.Modality,'SEG')
                    SEGheader = tmp;
                    SEG = squeeze(dicomread(elementFullFile));
                elseif strcmp(tmp.Modality,'MR') || strcmp(tmp.Modality,'PT') || strcmp(tmp.Modality,'CT')
                    sliceNumber = sliceNumber + 1;
                    volume{sliceNumber} = dicomread(elementFullFile); %double(dicomread(elementFullFile))
                    dicomHeaders = appendStruct(dicomHeaders,tmp); % appendStruct()
                    % refer <https://github.com/mvallieres/radiomics/>.
                end
            end
        end
        %    if waitB
        %        waitbar(elementNumber/nElements,waitbarHandle);
        %    end
    end % end reading for one patient
    nSlices = sliceNumber; % Total number of slices
    volume = volume(1:nSlices); % Suppress empty cells in images
    
    % DETERMINE THE SCAN ORIENTATION
    dist = [abs(dicomHeaders(2).ImagePositionPatient(1) - dicomHeaders(1).ImagePositionPatient(1)), ...
        abs(dicomHeaders(2).ImagePositionPatient(2) - dicomHeaders(1).ImagePositionPatient(2)), ...
        abs(dicomHeaders(2).ImagePositionPatient(3) - dicomHeaders(1).ImagePositionPatient(3))];
    [~,index] = max(dist);
    if index == 1
        orientation = 'Sagittal';
    elseif index == 2
        orientation = 'Coronal';
    else
        orientation = 'Axial';
    end
    
    % For original DICOM image series, SORT THE IMAGES AND DICOM HEADERS after
    % extracting the right order indices
    slicePositions = zeros(1,nSlices);
    for sliceNumber = 1:nSlices
        slicePositions(sliceNumber) = dicomHeaders(sliceNumber).ImagePositionPatient(index);
    end
    %[~,indices_img] = sort(slicePositions, 'descend'); % default is ascending, will
    % get indices_img like 8,7,6,5,4,3,2,1, inverse to the DICOM files we have and
    % might be inconsist with other tasks like dicom_read_volume()
    [~,indices_img] = sort(slicePositions);
    
    volume = cell2mat(volume(indices_img));
    dicomHeaders = dicomHeaders(indices_img);
    
    % transform values stored as uint16 to int16(Housefield Units, HU)
    type = dicomHeaders(1).Modality;
    if strcmp(type,'PT') || strcmp(type,'CT')
        if strcmp(type,'PT')
            type = 'PET';
        end
        for i=1:size(volume,3)
            volume(:,:,i)=volume(:,:,i)*dicomHeaders(i).RescaleSlope + dicomHeaders(i).RescaleIntercept;
        end
    end
    
    if exist('SEGheader', 'var') && isstruct(SEGheader)
        % For DSO derived mask volume, SORT THE IMAGES AND DICOM HEADERS after
        % extracting the right order indices
        % get and extract image positions
        slicePositions_DSO = zeros(1,nSlices);
        for nSDSO = 1:nSlices
            tmpSDSO = SEGheader.ReferencedSeriesSequence. ...
                Item_1.ReferencedInstanceSequence.(['Item_' num2str(nSDSO)]). ...
                ReferencedSOPInstanceUID;
            
            slice_index = structfind(dicomHeaders,'SOPInstanceUID',tmpSDSO); % find the slice index in dicomHeaders to identify the img slice corresponding to this DSO slice
            slicePositions_DSO(nSDSO) = dicomHeaders(slice_index).ImagePositionPatient(index);
        end
        %[~,indices_DSO] = sort(slicePositions_DSO, 'descend'); % default is ascending, will
        % get indices_img like 8,7,6,5,4,3,2,1, inverse to the DICOM files we have
        [~,indices_DSO] = sort(slicePositions_DSO);
        
        % NOTE: slicePositions and slicePositions_DSO should have the same values
        % set, but probably have different order of the values. That's why we need
        % to sort them separately to get the dicom images volume and the mask
        % volume with the same slice orders.
        
        SEGvol = SEG(:,:,indices_DSO);
        
        dicomInfo(j).SEGvol = SEGvol;
        dicomInfo(j).SEGheader = SEGheader;
    end
    
    % calculate slice spacing
    s1 = round(0.5*nSlices); s2 = round(0.5*nSlices) + 1; % Slices selected to
    % calculate slice spacing
    SliceSpacing =  sqrt(sum((dicomHeaders(s1).ImagePositionPatient - dicomHeaders(s2).ImagePositionPatient).^2));
    dicomInfo(j).volume = volume;
    
    dicomInfo(j).headers = dicomHeaders;
    dicomInfo(j).dicomInfo1slice = dicomHeaders(1);
    dicomInfo(j).SliceSpacing = SliceSpacing;
    dicomInfo(j).SliceThickness = dicomHeaders(1).SliceThickness;
    if isfield(dicomHeaders,'RescaleType')
        dicomInfo(j).RescaleType = dicomHeaders(1).RescaleType; % not all DICOM
    % series has this fieldname. TCGA collection "TCGA-HNSC" does, while
    % "NSCLC-Radiogenomics" does not.
    end
    dicomInfo(j).RescaleSlope = dicomHeaders(1).RescaleSlope;
    dicomInfo(j).RescaleIntercept = dicomHeaders(1).RescaleIntercept;
end % end reading all patients
% dicomData will contain data from 2 patients.

%% check if the patIDs are unique.
patIDs = {dicomInfo.patID};
if length(unique(patIDs)) == size(dicomInfo,2)
    fprintf('All patIDs are unique, ready to move on\n');
else fprintf('One or more patIDs are duplicated, PLEASE correct the ID\n');
end;



%% start to extract radiomic features
for i = 1:length(All_Data_Path) % 
    fprintf('Start to process %d th case \n',i);
    
    % to make sure DataVolume and MaskVolume(derived from DSO file) have the same slices orders,
    % recalculate DataVolume based on 
    DataVolume = dicomInfo(i).volume; 
    % use braces {} to extract the contents instead of a cell from
    % "dicomData" with parentheses ()
    
    
    if isfield(dicomInfo, 'SEGheader') && isstruct(dicomInfo(i).SEGheader) %exist('dicomInfo(i).SEGheader','var'), i is dynamic and is not accepted in exist().
        MaskVolume = dicomInfo(i).SEGvol;
    else
        MaskVolume = DataVolume;
        MaskVolume(MaskVolume~=0)=1;
    end
    % create a mask to pass to Func_RadiomicFeatureExtraction_Pipeline_v2() 
    % , a 2D or 3D array of dimensions corresponding to 'volume'. The mask 
    % contains 1's in the region of interest (ROI), and 0's elsewhere.
    
%     do_visualization=0;    
%     if do_visualization
%         save_folder = 'temp_file/';
%         name_case = i;
%         SavePath = [save_folder strcat(int2str(name_case), '_Headneck.avi')];
%         Func_Tumor_Visualization(DataVolume,MaskVolume,SavePath); % Func_Tumor_Visualization() not found. By Chao.
%     end
    
    %feature extraction
    %--------------------------------------------------------------------------
    % These parameters are mostly for prepareVolume(volume,mask,scanType,
    % pixelW,sliceS,R,scale,textType,quantAlgo,Ng) from 
    % <https://github.com/mvallieres/radiomics/>
    % The values assigments were modified by Chao according to the original
    % DICOM info embedded in the DICOM files.
    Para.type = 'CT'; % scanType: equals to prepareVolume()'s 'scanType': String 
    % specifying the type of scan analyzed. Either 'PETscan', 'MRscan' or 
    % 'Other'.
    Para.pixel_space = dicomInfo(i).headers(1).PixelSpacing(1);% pixelW: pixel width/spacing, Numerical value specifying the 
    % in-plane resolution (mm) of 'volume'.
    Para.slice_space = dicomInfo(i).SliceSpacing; % sliceS: Numerical value specifying the slice spacing (mm) of 'volume'.
    % Slice Thickness differs from Slice Spacing: more see https://stackoverflow.com/questions/14930222/how-to-calculate-space-between-dicom-slices-for-mpr
    Para.R_mat = 1;% R: Numerical value specifying the ratio of weight to 
    % band-pass coefficients over the weigth of the rest of coefficients 
    %(HHH and LLL). Provide R=1 to not perform wavelet band-pass filtering. 
    Para.scale_cell = 1;% scale: Numerical value specifying the scale at which 
    % 'volume' is isotropically resampled (mm). If a string 'pixelW' is 
    % entered as input, the volume will be isotropically resampled at the 
    % initial in-plane resolution of 'volume' specified by 'pixelW'. If 
    % scale_cell = 1, is isotropically resampled to a voxel size of 1x1x1
    % mm3.
    Para.textType = ''; % textType: String specifying for which type of textures 
    % 'volume' is being prepared. Either 'Global' or 'Matrix'. If 'Global',
    % the volume will be prepared for Global texture features computation. 
    % If 'Matrix',the volume will be prepared for matrix-based texture 
    % features computation (i.e. GLCM, GLRLM, GLSZM, NGTDM).
    % NOTE: textType will not be used later. textType will specifically assigned
    % with a string as 'Global' or 'Matrix' to prepareVolume() if
    % necessary.
    Para.algo_cell = 'Equal';% rescale into 0 - 255?. quantAlgo: String specifying the
    % quantization algorithm to use on 'volume'. Either 'Equal' for 
    % equal-probability quantization, 'Lloyd' for Lloyd-Max quantization, 
    % or 'Uniform' for uniform quantization. Use only if textType is set to
    % 'Matrix'.
    Para.Ng_mat = 256; % Ng: Integer specifying the number of gray levels in 
    % the quantization process. Use only if textType is set to 'Matrix'.
   
    fprintf('Extracting features for %d th case...\n', i);
    
    radiomicFeatures_tmp = radiomicFeatureExtraction_Pipeline(DataVolume,MaskVolume,Para);
    radiomicFeatures_tmp.patID = dicomInfo(i).headers.PatientID;
    if i == 1
        radiomicFeatures = radiomicFeatures_tmp;
    else
        radiomicFeatures = [radiomicFeatures,radiomicFeatures_tmp];
    end
    
    fprintf('Done the %d th case \n\n\n',i);
    DataVolume=[];
    MaskVolume=[];    
end
% end

%% Save radiomic features % modified by Chao
SaveFilePath = [outputDir prefix_of_output_csv datestr(now,formatOut) '.mat'];
save(SaveFilePath,'radiomicFeatures');

% change struct fieldnames, i.e. make the radiomic features structured as
% groupname_featureName. The grouping rules are based on the supplement of Chao's
% HNSC_Molecular_types paper.
oldNames = fieldnames(radiomicFeatures);
for nameIdx = 1:numel(oldNames)
    if ~strcmp(oldNames{nameIdx},'patID')
        if startsWith(oldNames{nameIdx}, 'firstOrder')
            newName = ['intensity_based_General_' extractAfter(oldNames{nameIdx},10)];
        elseif startsWith(oldNames{nameIdx}, 'GaborBank')
            newName = ['filter_based_Gabor_' extractAfter(oldNames{nameIdx},9)];
        elseif startsWith(oldNames{nameIdx}, 'GLCM')
            newName = ['textural_GLCM_' extractAfter(oldNames{nameIdx},4)];
        elseif startsWith(oldNames{nameIdx}, 'GLRLM')
            newName = ['textural_GLRLM_' extractAfter(oldNames{nameIdx},5)];
        elseif startsWith(oldNames{nameIdx}, 'GLSZM')
            newName = ['textural_GLSZM_' extractAfter(oldNames{nameIdx},5)];
        elseif startsWith(oldNames{nameIdx}, 'HistBin')
            newName = ['intensity_based_Histogram_' extractAfter(oldNames{nameIdx},7)];
        elseif startsWith(oldNames{nameIdx}, 'HOGvalue')
            newName = ['textural_HOG_value_' extractAfter(oldNames{nameIdx},8)];
        elseif startsWith(oldNames{nameIdx}, 'LBPvalue')
            newName = ['textural_LBP_value_' extractAfter(oldNames{nameIdx},8)];
        elseif startsWith(oldNames{nameIdx}, 'NGTDM')
            newName = ['textural_NGTDM_' extractAfter(oldNames{nameIdx},5)];
        elseif startsWith(oldNames{nameIdx}, 'shapeSize')
            newName = ['shape_size_' extractAfter(oldNames{nameIdx},9)];
        elseif startsWith(oldNames{nameIdx}, 'wavelet')
            newName = ['filter_based_Wavelet_' extractAfter(oldNames{nameIdx},7)];
        end
        
        [radiomicFeatures.(newName)] = radiomicFeatures.(oldNames{nameIdx});
        radiomicFeatures = rmfield(radiomicFeatures,oldNames{nameIdx});
    end
end

% transform struct to table
radiomicFeaturesTable = struct2table(radiomicFeatures);

% move the last column (patID) to the first as conventional practice
% NO-NEED: radiomicFeaturesTable = [radiomicFeaturesTable(:,size(radiomicFeaturesTable,2)) radiomicFeaturesTable(:, 1:(size(radiomicFeaturesTable,2)-1))];

SaveTablePath = [outputDir prefix_of_output_csv datestr(now,formatOut) '.csv'];
writetable(radiomicFeaturesTable, SaveTablePath);

features_extracted = radiomicFeaturesTable;
end