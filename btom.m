% btom - takes output from cellPACK and turns it into a synthetic tomogram .mrc file

tic; clear %bbchange

name = 'output'; % INPUT - locations and rotations of all particles. CHANGE as needed

[~, work_dir, ~] = fileparts(pwd);
if ~strcmp(work_dir,'tomosimu_sandbox')
    cd ..
end
cd code; work_dir=pwd;

cellPACK=importdata(strcat(name,'.txt'));

protein_names=unique(cellPACK.textdata);
%protein_numbers=cellPACK.data(:,13);
locations=cellPACK.data(:,1:3);

rotations=cell(numel(cellPACK.data(:,1)),1);
for i=1:numel(cellPACK.data(:,1))
    rotations{i}=[[cellPACK.data(i,4), cellPACK.data(i,5), cellPACK.data(i,6)]; [cellPACK.data(i,7), cellPACK.data(i,8), cellPACK.data(i,9)]; [cellPACK.data(i,10), cellPACK.data(i,11), cellPACK.data(i,12)]];
end

bounding_box=[ceil(max(locations(:,1)))+20 ceil(max(locations(:,2)))+20 ceil(max(locations(:,3)))+20];
%bounding_box=[100 100 100]; %CHANGE to input bounding box from file!!!!

tomogram_size = bounding_box; % in pixels? nm? I think this is just the number of boxes in the tomogram.

disp('Initialization Start');
ws = struct();
ws.map.map_resolution_s = '40'; %bb This (cubed) is the size of all of the individual protein maps. They will all be the same size. Situs generates maps of different sizes, so they need to be rescaled
ws.map.map_resolution = str2double(ws.map.map_resolution_s);
ws.map.rescale_individual_maps = false;
ws.map.protein_names = protein_names;
ws.sys.data_dir = strcat(work_dir,'/PDB');
ws.map.volumes =  GenerateSimulationMap.load_maps(ws); %density maps for pdbs

% resize so all same - why?
[vols_large, max_siz] = VolumeUtil.resize_vols(ws.map.volumes, 1.0);

% definition of reconstruction model parameters
% set values for SNR, missing wedge and ctf
ws.reconstruction_param.model = struct();
ws.reconstruction_param.model.SNR = 1; % was 0.05. 
ws.reconstruction_param.model.missing_wedge_angle = 30;
ws.reconstruction_param.model.ctf = GenerateSimulationMap.get_ctf_param(ws.map.map_resolution);
ws.reconstruction_param.model.ctf.voltage=300;

% set the contour level for each macromolecular complex
contour_threshold = 0.2;
for i=1:numel(vols_large)
        density_threshold(i)=max(vols_large{i}(:))*contour_threshold;
        standard_sphere(i)=amprs_wthr(vols_large{i},density_threshold(i),max_siz(1)); %GET RID OF!!! USE CENTROIDS FROM CELLPACK? THIS IS WHERE THE VARIABILITY COMES FROM. IN amprs_whtr IS FAULTY minboundsphere. Actually, may be OK - only problem was when using radii to identify proteins. using for rotation centers may still be ok. But there might be a more efficient way to do it.
        %complex_volume(i)=sum(vols_large{i}(:)>density_threshold(i)); %USELESS?
end

% add each particle to its corresponding location to generate simulated tomogram
tic
vol_den=zeros(tomogram_size);
disp('adding particle            ')
for i=1:size(locations,1)
    protein_number=find(strcmp(cellPACK.textdata(i),protein_names));
    fprintf(1,'\b\b\b\b\b\b\b\b\b\b%10.0f',i);
    rotation_center = standard_sphere(protein_number).c;
    shifting_to_location = locations(i,1:3);
    vol_t_mut=VolumeUtil.rotate_vol_pad0(vols_large{protein_number}, rotations{i}, rotation_center,shifting_to_location,tomogram_size,'cubic');
    vol_den=vol_den+vol_t_mut;
end
fprintf('\n'); toc

% apply back projection to realistically simulate tomogram
vol_den_bp=GenerateSimulationMap.backprojection_reconstruction(ws.reconstruction_param, vol_den, ws.reconstruction_param.model.SNR);

cd ../tomograms
tom_mrcwrite(vol_den_bp,'name',strcat(name,'.mrc'),'style','fei');
cd ..; toc

