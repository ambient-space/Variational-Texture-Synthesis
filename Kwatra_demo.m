clear all
addpath(genpath('./utility/'));
addpath('./../arbyreed textures/');

figure
files = dir('./../arbyreed textures/*.png');

batch_name = './../test/';
mkdir(batch_name);


first_file = 'ghost_peppers.png';
b = 1;
for file = files'
    
    if b && ~strcmp(file.name,first_file)
        continue;
    else
        b =0;
    end
    
params.in_scale = 1;
params.out_scale = 1;
params.iter = 12;

disp(file.name)
x0 = load_texture(file.name);
[m,n,~] = size(x0);
if m == n
    x0 =imresize(x0,[256,256]);
else
x0 =imresize(x0,[480,640]);
end
x0 = Spectrum.periodic(x0);

params.scales  = [1,.5,.25];
t1 = Texture(x0,params);

wc_params.num_levels = 4;
t1.add_constraint(WaveletConstraint(t1,wc_params));
% t1.add_constraint(HistogramConstraint(t1));

mc_params.blocksize = 3;
mc_params.x_dataratio = .4;
mc_params.y_dataratio = .5;
mc_params.num_levels = 4;
mc_params.scales = params.scales;

% t1.add_constraint(MRFConstraint(t1,mc_params));
t1.add_constraint(WaveletAncestryMRFConstraint(t1,mc_params));


t1.run_variational_synthesis()


info = Analyze.get_tiling(t1.x0,t1.y, 7 , 0);

synth_folder = [batch_name,'synthesis/'];
tiling_folder = [batch_name,'tilings/'];

mkdir(synth_folder);
mkdir(tiling_folder);

imwrite(t1.y,[synth_folder,file.name(1:end-4),'.png'],'png');
imwrite(info.tiling_image,[tiling_folder,file.name(1:end-4),'.png'],'png');

if exist([batch_name,'innovation_capacities.mat'])
    innovation_capacities = load([batch_name,'innovation_capacities.mat']);
end

% add innovation capacity to struct, then save
eval(['innovation_capacities.',file.name(1:end-4),' = ',num2str(info.innovation_capacity)]);


    
save([batch_name,'innovation_capacities.mat'],'-struct','innovation_capacities');

end