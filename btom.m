% btom - takes output from cellPACK and turns it into a synthetic tomogram .mrc file
% "loadjson" requires "JSONlab_ a toolbox to encode_decode JSON files"
% Toolbox
% "RotationMatrix" and "quaternion" require "quaternion" Toolbox
% Also dependent on TOM_Release_2008 Toolbox
% Also dependent on tomosimu software from Alber lab

tic; clear

name = 'RECIPE-Alber_model_ALL-1QO1-2REC_500_res_tr.json'; % INPUT - cellPACK result _tr .json file

[~, work_dir, ~] = fileparts(pwd);
if ~strcmp(work_dir,'tomosimu_sandbox')
    cd ..
end
cd code; 
work_dir=pwd;

recipe = loadjson(name);
compartments = fieldnames(recipe);
%ingredients=fieldnames(recipe.(compartments{2})); %NEED TO EXPAND for multiple compartments
%protein_names=unique(cellPACK.textdata);
%protein_numbers=cellPACK.data(:,13);
%locations=cellPACK.data(:,1:3);

%bounding_box=[ceil(max(locations(:,1)))+20 ceil(max(locations(:,2)))+20 ceil(max(locations(:,3)))+20];
bounding_box=[500 500 200]; %CHANGE to input bounding box from file!!!!

tomogram_size = bounding_box; % in pixels? nm? I think this is just the number of boxes in the tomogram.

for i=1:numel(compartments)-1
    ingredients=fieldnames(recipe.(compartments{i+1}).ingredients);
%    pdbs=cell(numel(ingredients),1);
    for j=1:size(ingredients)
        pdbs{i,j}=recipe.cytoplasme.ingredients.(ingredients{j}).source.pdb;
    end
end
protein_names=unique(pdbs);

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
ws.reconstruction_param.model.SNR = 50; % was 0.05. 
ws.reconstruction_param.model.missing_wedge_angle = 30;
ws.reconstruction_param.model.ctf = GenerateSimulationMap.get_ctf_param(ws.map.map_resolution);
ws.reconstruction_param.model.ctf.voltage=300;

% add each particle to its corresponding location to generate simulated tomogram
rotation_center = ws.map.map_resolution/2 * [1 1 1]; %MUST CHANGE - rotates around center of .situs map. cellPACK rotates around center of mass (sort of - actually average positions of atoms, not weighted). They must be made to agree. Calculate center of mass here?
tic
vol_den=zeros(tomogram_size);

disp(strcat('recipe',recipe.recipe.name))

for i=2:numel(compartments)
    disp(strcat('compartment=',compartments{i}))
    ingredients=fieldnames(recipe.(compartments{i}).ingredients);
    for j=1:numel(ingredients)
        disp(strcat('ingredient=',ingredients{j}))
        pdb=recipe.cytoplasme.ingredients.(ingredients{j}).source.pdb;
        disp(strcat('pdb=',pdb))
        particles = recipe.cytoplasme.ingredients.(ingredients{j}).results;
        disp('adding particle            ')
        for k=1:numel(particles)
            protein_number=find(strcmp(pdb,protein_names));
            fprintf(1,'\b\b\b\b\b\b\b\b\b\b%10.0f',k);
            shifting_to_location = particles{k}{1}/10;
            rotation_matrix=RotationMatrix(quaternion(particles{k}{2}));
            vol_t_mut=VolumeUtil.rotate_vol_pad0(vols_large{protein_number}, rotation_matrix, rotation_center,shifting_to_location,tomogram_size,'cubic');  %this is the time-taker!
            vol_den=vol_den+vol_t_mut;
        end
    end
end
fprintf('\n'); toc

% apply back projection to realistically simulate tomogram
vol_den_bp=GenerateSimulationMap.backprojection_reconstruction(ws.reconstruction_param, vol_den, ws.reconstruction_param.model.SNR);
percent_volume_filled = 100 - sum(vol_den(:)==0)/numel(vol_den);
disp(strcat('percent volume filled=',percent_volume_filled))

cd ../tomograms
tom_mrcwrite(vol_den_bp,'name',strcat(name,'.mrc'),'style','fei');
cd ..; toc

