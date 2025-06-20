# HR-pQCT-Erosion-Detection    
The Matlab implementation of "Automatic 3D joint erosion detection for the diagnosis and monitoring of rheumatoid arthritis using hand HR-pQCT images", CMIG, 2023.   
The test of this codes needs the HR-pQCT joint volume. For the enquiry of input joint volume, you can contact Isaac Cheng (ichength@cuhk.edu.hk) and Xuechen Zhang (xzhangem@connect.ust.hk).


### Descriptions 
The codes are tested on MATLAB R2019 or above version.

Cropping tool takes DICOM slices of whole hand joints, and outputs the individual 3D metacarpophalangeal (MCP) in .nii format. 

The levelset_utils file contains the codes for level set-based cortical surface construction, and the curvature_utils file contains the codes for surface curvature-based features for erosion detection. UI.mlapp is a simple MATLAB UI with input of MCP volume (in .nii format) and outputs the bone segmentation and erosion detection results (in .nii format).    
