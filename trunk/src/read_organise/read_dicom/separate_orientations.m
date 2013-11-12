function ACS = separate_orientations(patient)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%  Separate the different views: axial, coronal and sagittal for a given
%%  patient
%%
%%  Inputs:  - patient -> struct containing 3 cells (images, dicom_inf, orientation)
%%  Outputs: - ACS -> struct containing 3 matrices (axial, coronal, sagittal)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cont_a = 1;
cont_c = 1;
cont_s = 1;
for i=1:size(patient.images,2)
    
    if(strcmp(patient.orientation{i},'A'))
        ACS.axial(:,:,cont_a) = patient.images{i};
        cont_a = cont_a + 1;
    elseif(strcmp(patient.orientation{i},'C'))
        ACS.coronal(:,:,cont_c) = patient.images{i};
        cont_c = cont_c + 1;
    elseif(strcmp(patient.orientation{i},'S'))
        ACS.sagittal(:,:,cont_s) = patient.images{i};
        cont_s = cont_s + 1;
    end
    
end