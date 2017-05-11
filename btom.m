% btom - takes output from cellPACK and turns it into a synthetic tomogram .mrc file
% "loadjson" requires "JSONlab_ a toolbox to encode_decode JSON files" Toolbox
% "RotationMatrix" and "quaternion" require "quaternion" Toolbox by Mark
% Tincknell - https://sourceforge.net/projects/qtfm/  qtfm_2_5.zip
% Also dependent on TOM_Release_2008 Toolbox
% Also dependent on tomosimu software from Alber lab

tic; clear

model_dir = '/Users/mac/Documents/Alber_model_3/';
model_name = 'RECIPE-Alber_model_ALL_res_tr.json'; % INPUT - cellPACK result _tr .json file
model_path = [model_dir model_name];
out_name = [model_name ''];
%bounding_box=[ceil(max(locations(:,1)))+20 ceil(max(locations(:,2)))+20 ceil(max(locations(:,3)))+20];
bounding_box=[100 100 100]; %CHANGE to input bounding box from file!!!!

% make sure we're in the right directory
[~, curr_dir, ~] = fileparts(pwd);
if ~strcmp(curr_dir,'btom')
    cd ..
end
cd code; 

recipe = loadjson(model_path);
compartments = fieldnames(recipe);

tomogram_size = bounding_box; % in pixels? nm? I think this is just the number of boxes in the tomogram.

for i=2:numel(compartments)
    ingredients=fieldnames(recipe.(compartments{i}).ingredients);
%    pdbs=cell(numel(ingredients),1);
    for j=1:size(ingredients)
        pdbs{i-1,j}=recipe.cytoplasme.ingredients.(ingredients{j}).source.pdb;
    end
end
pdbs=unique(pdbs); % this will remove duplicates and rearrange the pdbs into alphabetical order


cd ..


%make situs files
pdb2vol_loc = './Situs_2.8/src/pdb2vol';
for i=1:numel(pdbs)
    pdbfile = [model_dir,'PDB/',pdbs{i},'/',pdbs{i},'.pdb'];
    situsfile = [model_dir, 'PDB/',pdbs{i},'/40','.situs'];
    settingsfile = './Situs_2.8/settings.txt';
    system_string = [pdb2vol_loc ' ' pdbfile ' ' situsfile ' < ' settingsfile];
    system(system_string)
end


disp('Initialization Start');
ws = struct();
ws.map.map_resolution_s = '40'; % This is the resolution of the maps generated by situs pdb2vol
ws.map.map_resolution = str2double(ws.map.map_resolution_s);
ws.map.rescale_individual_maps = false;
ws.map.protein_names = pdbs;
ws.sys.data_dir = strcat(model_dir,'PDB');
ws.map.volumes = GenerateSimulationMap.load_maps(ws); %density maps for pdbs

% resize so all same - why?
[vols_large, max_siz] = VolumeUtil.resize_vols(ws.map.volumes, 1.0);
%vols_large = ws.map.volumes;

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
        pdb=recipe.(compartments{i}).ingredients.(ingredients{j}).source.pdb;
        disp(strcat('pdb=',pdb))
        particles = recipe.(compartments{i}).ingredients.(ingredients{j}).results;
        disp('adding particle            ')
        protein_number=find(strcmp(pdb,pdbs));
        for k=1:numel(particles)
            fprintf(1,'\b\b\b\b\b\b\b\b\b\b%10.0f',k);
            shifting_to_location = particles{k}{1}/10;
            rotation_matrix=RotationMatrix(quaternion(particles{k}{2}));
            vol_t_mut=VolumeUtil.rotate_vol_pad0(vols_large{protein_number}, rotation_matrix, rotation_center,shifting_to_location,tomogram_size,'cubic');  %this is the time-taker!
            vol_den = vol_den + vol_t_mut;
        end
        fprintf('\n')
    end
end
fprintf('\n'); toc

% apply back projection to realistically simulate tomogram
vol_den_bp=GenerateSimulationMap.backprojection_reconstruction(ws.reconstruction_param, vol_den, ws.reconstruction_param.model.SNR);
percent_volume_filled = 100 - 100*sum(vol_den(:)==0)/numel(vol_den);
fprintf(1,'percent volume filled=%f\n',percent_volume_filled)

%cd ../tomograms
tom_mrcwrite(vol_den_bp,'name',[model_dir out_name '.mrc'],'style','fei');
toc
