# Colorectal_cancer_radiomics
Contains code to predict near-complete response to Chemoradiotherapy in colorectal cancers using radiomics on MR images.

The codebase is divided in two jupyter notebooks: one for preprocessing the data (Preprocessing with image alignment and b-value interpolation.ipnyb) and another for building predictive models (Model building.ipnyb). 

The preprocessing notebook handles the different preprocessing that was required for this work, including constructing DWI images for b-value 800 using linear interpolation on 3-5 volumes with different available b-values, creating ADC maps by combining several DWI volumes with different b-values and transferring the tumor ROIs from the T2-weighted images to corresponding DWI images and ADC maps.

The model building notebook was used for experimenting with different models based on radiomic features, choosing the best performing models on the training set and finally finetuning models and validating them on the validation set.

While both of the above steps were performed in Python, the intermediate step of actually extracting radiomics features from the MR volumes and the ROIs was done in MATLAB. Our in-house feature extractor is included here as a zipped file, which needs to be uncompressed before using. The entry-point to the feature extractor is through the matlab code file (mainFunc_radiomicFeatureExtraction.m). 
