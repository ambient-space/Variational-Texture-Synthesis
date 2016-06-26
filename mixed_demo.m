clear all
addpath(genpath('./utility/'));
addpath('./../arbyreed textures/');

figure
files = dir('./../arbyreed textures/*.png');

batch_name = './../arbyreed wa mrf batch test 3/';
mkdir(batch_name);

file1 = 'packaged_candy.png';
file2 = 'bolts.png';

b = 1; % 1 = use first file

    
params.out_scale = 1;
params.iter = 32;

disp(file1)
disp(file2)
x1 = load_texture(file1);
x1 =imresize(x1,[480,640]/2);
x2 = load_texture(file2);
x2 =imresize(x2,[480,640]/2);

x1 = Spectrum.periodic(x1);
x2 = Spectrum.periodic(x2);

params.scales  = [1,.5,.25];
params.band_weights = [.6,.3,.6,.6];
t1 = TextureMixture(x1,x2,params);

wc_params.num_levels = 4;
t1.add_constraint(MixedWaveletConstraint(t1,wc_params));
% t1.add_constraint(HistogramConstraint(t1));

mc_params.blocksize = 2;
mc_params.x_dataratio = .6;
mc_params.y_dataratio = .6;
mc_params.num_levels = 4;
% t1.add_constraint(SpectrumConstraint(t1));
% t1.add_constraint(MRFConstraint(t1,mc_params));

t1.add_constraint(MixedWaveletAncestryMRFConstraint(t1,mc_params));


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
