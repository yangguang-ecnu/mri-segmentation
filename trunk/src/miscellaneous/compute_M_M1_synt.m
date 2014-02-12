function [M,M_1,transform] = compute_M_M1_synt(patient, ortho, direction, ImagePositionPatient)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%  Given the DICOM information of the patient, compute the transform from
%%  the coordinate frame to the reference coordinate system (RCS)
%%
%%  Inputs:  1. patient -> a struct containing the DICOM of an image
%%           2. ortho   -> 0/1 assume the directions are orthonormal/ the DICOM info
%%  Outputs: 1. M       -> 4x3 matrix
%%           2. M_1     -> 3x4 matrix, the inverse of M
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

transform = eye(3);

spacing = double((patient.PixelSpacing(1) * patient.SpacingBetweenSlices * patient.ImagesInAcquisition))/double((patient.Rows*patient.PixelSpacing(1)));%0.3094;;;

tmp_A = zeros(4,3);
tmp_A(4,3) = 1;

if ortho
    if direction == 1
        
        tmp_A(1:3,1) = [1; 0; 0].*patient.PixelSpacing(1);
        tmp_A(1:3,2) = [0; 1; 0].*patient.PixelSpacing(2);
        tmp_A(3,3)   = 1;
        
        cross_p   = cross(patient.ImageOrientationPatient(1:3),patient.ImageOrientationPatient(4:6));
        r         = vrrotvec(cross_p, [0 0 1]');
        transform = vrrotvec2mat(r);
        
        new_origin = transform * patient.ImagePositionPatient;


    end
    if direction == 2
        tmp_A(1:3,1) = [0; 1;  0] .* patient.PixelSpacing(2);
        tmp_A(1:3,2) = [0; 0; -1] .* spacing;%patient.SpacingBetweenSlices;
        tmp_A(1,3)   = 1;
        
%         cross_p = cross(patient.ImageOrientationPatient(1:3),patient.ImageOrientationPatient(4:6));
%         r = vrrotvec(cross_p, [-1 0 0]');
%         transform = vrrotvec2mat(r);
        
%         cross_p   = cross(patient.ImageOrientationPatient(1:3),patient.ImageOrientationPatient(4:6));
%         r         = vrrotvec(cross_p, [0 0 1]');
%         transform = vrrotvec2mat(r);
        
        new_origin = transform * ImagePositionPatient;
    
    end
    if direction == 3
        tmp_A(1:3,1) = [-1; 0;  0] .* patient.PixelSpacing(1);
        tmp_A(1:3,2) = [ 0; 0; -1] .* spacing;%
        tmp_A(2,3) = 1;
        
%         cross_p = cross(patient.ImageOrientationPatient(1:3),patient.ImageOrientationPatient(4:6));
%         r = vrrotvec(cross_p, [0 1 0]');
%         transform = vrrotvec2mat(r);
        
%         cross_p   = cross(patient.ImageOrientationPatient(1:3),patient.ImageOrientationPatient(4:6));
%         r         = vrrotvec(cross_p, [0 0 1]');
%         transform = vrrotvec2mat(r);
        
        new_origin = transform * ImagePositionPatient;
    end
else
    tmp_A(1:3,1) = patient.ImageOrientationPatient(1:3).*patient.PixelSpacing(1);
    tmp_A(1:3,2) = patient.ImageOrientationPatient(4:6).*patient.PixelSpacing(2);
    
    new_origin = transform * patient.ImagePositionPatient;
end

%%
% patient.ImagePositionPatient
% new_origin = transform * patient.ImagePositionPatient
tmp_b = eye(4);
tmp_b(1:3,4) = new_origin;
% tmp_b(1:3,4) = patient.ImagePositionPatient;


% M = [transform zeros(3,1);0 0 0 1] * tmp_b * tmp_A
% M_1 = pinv(tmp_A) * pinv(tmp_b) * pinv([transform zeros(3,1);0 0 0 1]);

M =  tmp_b * tmp_A;
M_1 = pinv(tmp_A) * pinv(tmp_b);


% tmp_A
% tmp_b