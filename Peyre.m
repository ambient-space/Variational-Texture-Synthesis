clear all
close all
addpath(genpath('./utility/'));
addpath('./../arbyreed textures/');

figure
files = dir('./../arbyreed textures/*.png');

batch_name = './../arbyreed peyre batch 1/';
mkdir(batch_name);


first_file = 'colored_animals.png';
b = 1; % 1 = use first file
for file = files'
    
    if b && ~strcmp(file.name,first_file)
        continue;
    else
        b =0;
    end
    
params.out_scale = 1;
params.iter = 32;

disp(file.name)
x0 = load_texture(file.name);
[m,n,~] = size(x0);
if m == n
    x0 =imresize(x0,[256,256]);
else
x0 =imresize(x0,[480,640]/2);
end
x0 = Spectrum.periodic(x0);

params.scales  = [1,.5,.25];
t1 = Texture(x0,params);

wc_params.num_levels = 4;

sc_params.blocksize = 10;
sc_params.x_dataratio = .3;
sc_params.y_dataratio = .25;
sc_params.sparsity = 4;
sc_params.dictsize = 2^10;

% t1.add_constraint(SparsityConstraint(t1,sc_params));
% t1.add_constraint(SparsityConstraint(t1,sc_params));
t1.add_constraint(SpectrumConstraint(t1));
t1.add_constraint(HistogramConstraint(t1));
t1.add_constraint(SparsityConstraint(t1,sc_params));

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