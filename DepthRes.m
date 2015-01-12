function DRes = DepthRes(SPC)
% This function extracts the location of the resolvable depth planes for
% a PC2 camera.
% It considers only the central row of the pixels on the image sensor.
% The intersection of the projected pixel with the start point of the light
% cones is extracted.
% In this implementation, the main lens in discarded.
% We look at the depth planes located after the focus of microlenses.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% camera parameters:

% f = 25;
% f_num = 1.8;
% pitch = 13.898614883422850808;
% n_x = 10;
% g = 25.1;
% acceptance_angle = 40; % Light acceptance angle of the image sensor pixels

% calculate the line equation for the proj_pix_size in the positive x (line p_projpix).
% sort the light cones based on their x.
% For all the positive x:
% 1-write the line equation for the start of the light cone (line p_LC),
% 2-find the intersection between line p_projpix and p_LC.
% 3-if the intersection point is valid (positive and not NaN) add it to DRes.

% Line p_projpix
% Pixel size on the sensor: pitch / n_x
% First point: lower extent of the central pixel: [x,z] = [-(pitch / n_x / 2), 0]
% Second point: Center of the lens array: [x,z] = [0, g]

% Line p_LC
% light cone number: i
% First point: light cone tip: [x,z] = [ALCm(i,4),ALCm(i,6)]
% Second point: [x,z] = [ ALCm(i,4) +  ALCm(i,7) ,  ALCm(i,6) + 1 ]

%fit linear polynomial

DRes = [];
load(SPC);

ALCm = LC_tempL_m_final;
ALCm = sortrows(ALCm,4);

p_projpix = polyfit( [0,g] , [ -( pitch / n_x / 2 ),0], 1);

for i = floor(size(ALCm,1)/2):size(ALCm,1)% floor(size(ALCm,1)/2)
    p_LC = polyfit([ALCm(i,6), ALCm(i,6) + 100] , [ALCm(i,4), ALCm(i,4) +  100*ALCm(i,7)] , 1);
    
    %calculate intersection
    z_intersect = fzero(@(z) polyval( p_projpix - p_LC, z ),3);
    %     x_intersect = polyval(p_projpix,z_intersect);
    
    if (z_intersect > 0) && (~isnan(z_intersect))
        DRes = [DRes, z_intersect];
    end
end

DRes = sort (DRes);
n_z = size(DRes);
zmin = min (DRes);
zmax = max (DRes);
figure;
plot (DRes,1,'*');
figure;
plot (DRes(1:(size(DRes,2)-1)),DRes(1:(size(DRes,2)-1))-DRes(2:size(DRes,2)))


DRname = ['DResPlanes_' num2str(n_l) 'x1_lenslet_output_decoupled_' num2str(n_x) 'pixel_p=' num2str(pitch/1000) 'mm_g=' num2str(g/1000) 'mm_f=' num2str(f/1000) 'mm_f#=' num2str(f_num) '.mat'];
save (DRname,'LC_tempL_m_final','pitch','n_l','g','n_x','f','f_num');

%--------------------------------------------------------------------------
% % Uncomment the following in the case of comparing the SPC results with 
% % the camera setup built in:
% % @inproceedings{langguth2014optical,
% %   title={Optical performance analysis of plenoptic camera systems},
% %   author={Langguth, Christin and Oberd{\"o}rster, Alexander and Br{\"u}ckner, Andreas and Wippermann, Frank and Br{\"a}uer, Andreas},
% %   booktitle={SPIE Optical Engineering+ Applications},
% %   pages={91920F--91920F},
% %   year={2014},
% %   organization={International Society for Optics and Photonics}
% % }
% %
% % The camera parameters in this setup are (from the above paper):
% %
% % f (microlenses) = 0.5mm
% % f-number (microlenses) = 4.8
% % F main lens = 16mm
% % Distance from the main lens to the plane of the SPC cone tips = 19.26mm 
% % microlens pitch = 0.195mm
% % microlens aperture 0.104mm
% % pixel size = 2.2um
% % number of pixels = 2592*1944
% % g = 0.4 + 0.125 + 0.3mm
% % number of pixels in a row behind one microlens = 89
% %
% % PC2 for this structure will be built with only 3 microlenses since the
% % paper wants to look at the depth resolution as the result of parallax 
% % between two adjacant microlenses.
%--------------------------------------------------------------------------
DRes_inrange = [];
for i=1:size(DRes,2)
    if DRes(1,i)>3960 && DRes(1,i)<5090
        DRes_inrange = [DRes_inrange , DRes(i)];
    end
end
s = 19260+1270+825+((19260+1270+825-DRes_inrange).*16000./((19260+1270+825-DRes_inrange)-16000));
for i = 1:size(s,2)
    S(1,i) = s(i);
    S(2,i) = i-1;
end
figure;
plot(S(2,:),S(1,:),'+');
% --------------------------------------------------------------------------
% Measurement results from the above paper
M(1,:) = 300000:50000:950000;
M(2,:) = [1.4   2.6     4   4.5     5.2     5.5     5.75    6.25    6.5     6.8     7.2     7.45    8.3     8.6];
M(3,:) = [0.8   0.5     0.5 0.4     0.5     0.45    0.55    0.7     0.5     0.25    0.2     0.2     0.35    0.5];
hold on
herrorbar(M(2,:),M(1,:),M(3,:),M(3,:));


