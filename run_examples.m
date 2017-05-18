
close all
clc

run('vlfeat-0.9.14/toolbox/vl_setup');

% Example 1: temple
img_n = 2; % The number of input images
in_name = cell(img_n,1);
in_name{1} = 'images/AANAP_01.png';
in_name{2} = 'images/AANAP_02.png';  

edge_list = [1,2]; % Each row represents an image pair to be aligned.
% In this example, there is only one pair (image 1 and image 2).
ref = 2; % the index of the reference image
tar = 1; % the index of the target image. Our warp is constructed from the homgraphy that maps from ref to tar
warp_type = 'ours'; % 'ours' for our warp. 'hom' for homography warp
zeroR_ON = 1; % whether we restrict the similarity warp in our warp to be no rotation (zeroR_ON=1) or not (zeroR_ON=0)
SPHP_stitching(in_name, edge_list, ref, tar, warp_type, zeroR_ON); % run stitching
% title('Example 1: temple. (Press any key to continue with the next example.)');
% pause
% close all
% 
% % Example 2: carpark
% img_n = 2;
% in_name = cell(img_n,1);
% in_name{1} = 'images/carpark_01.jpg';
% in_name{2} = 'images/carpark_02.jpg';
% edge_list = [1,2];
% SPHP_stitching(in_name, edge_list, 2, 1, 'ours', 1);
% title('Example 2: carpark. (Press any key to continue with the next example.)');
% pause
% close all
% 
% % Example 3: campus
% img_n = 2;
% in_name = cell(img_n,1);
% in_name{1} = 'images/campus_01.jpg';
% in_name{2} = 'images/campus_02.jpg';
% edge_list = [1,2];
% SPHP_stitching(in_name, edge_list, 2, 1, 'ours', 1);
% title('Example 3: campus. (Press any key to continue with the next example.)');
% pause
% close all
% 
% % Example 4: street
% img_n = 3;
% in_name = cell(img_n,1);
% in_name{1} = 'images/street_01.jpg';
% in_name{2} = 'images/street_02.jpg';
% in_name{3} = 'images/street_03.jpg';
% edge_list = [1,2; 2,3]; % Two image pairs. The 1st one is (I1, I2), the 2nd one is (I2, I3).
% SPHP_stitching(in_name, edge_list, 3, 1, 'ours', 1);
% title('Example 4: street. (Press any key to continue with the next example.)');
% pause
% close all
% 
% % Example 5: disneysea
% img_n = 2;
% in_name = cell(img_n,1);
% in_name{1} = 'images/disneysea_01.jpg';
% in_name{2} = 'images/disneysea_02.jpg';
% edge_list = [1,2];
% SPHP_stitching(in_name, edge_list, 2, 1, 'ours', 0);
% title('Example 5: disneysea. (Press any key to continue with the next example.)');
% pause
% close all
% 
% % Example 6: bridge
% img_n = 3;
% in_name = cell(img_n,1);
% in_name{1} = 'images/bridge_01.jpg';
% in_name{2} = 'images/bridge_02.jpg';
% in_name{3} = 'images/bridge_03.jpg';
% edge_list = [1,2; 2,3];
% SPHP_stitching(in_name, edge_list, 3, 1, 'ours', 0);
% title('Example 6: bridge. (Press any key to continue with the next example.)');
% pause
% close all
% 
% % Example 7: railtracks
% img_n = 2;
% in_name = cell(img_n,1);
% in_name{1} = 'images/railtracks_01.jpg';
% in_name{2} = 'images/railtracks_02.jpg';
% edge_list = [1,2];
% SPHP_stitching(in_name, edge_list, 2, 1, 'ours', 1);
% title('Example 7: railtracks. (Press any key to continue with the next example.)');
% pause
% close all
% 
% % Example 8: park
% img_n = 2;
% in_name = cell(img_n,1);
% in_name{1} = 'images/park_01.jpg';
% in_name{2} = 'images/park_02.jpg';
% edge_list = [1,2];
% SPHP_stitching(in_name, edge_list, 2, 1, 'ours', 1);
% title('Example 8: park. (Press any key to continue with the next example.)');
% pause
% close all
% 
% % Example 9: garden
% img_n = 2;
% in_name = cell(img_n,1);
% in_name{1} = 'images/garden_01.jpg';
% in_name{2} = 'images/garden_02.jpg';
% edge_list = [1,2];
% SPHP_stitching(in_name, edge_list, 2, 1, 'ours', 1);
% title('Example 9: garden. (Press any key to continue with the next example.)');
% pause
% close all
% 
% % Example 10: building
% img_n = 3;
% in_name = cell(img_n,1);
% in_name{1} = 'images/building_01.jpg';
% in_name{2} = 'images/building_02.jpg';
% in_name{3} = 'images/building_03.jpg';
% edge_list = [1,2; 2,3];
% SPHP_stitching(in_name, edge_list, 3, 1, 'ours', 1);
% title('Example 10: building. (This is the last example.)');
