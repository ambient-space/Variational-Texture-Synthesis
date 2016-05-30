clear all
addpath('./../texture synthesis/textures/');
addpath(genpath('./utility/'));


rng(7172)

figure
params.in_scale = 1;
params.out_scale = 2;
params.iter = 8;
% x0 = load_texture('wood.png');
% x0 = load_texture('flowers_vis.png');
% x0 = load_texture('raddish_simon.jpg');
%  x0 = load_texture('shells_cc.jpg');
% x0 = load_texture('canvas_cc.jpg');
% x0 = load_texture('pebbles_cc.jpg');
% x0 = load_texture('marble_cg.jpg');
% x0 = load_texture('frozen_cc.jpg');
x0 = load_texture('lava5_cc.jpg');
% x0 = load_texture('sweet_potato_cc.jpg');

x0 =imresize(x0,[256,256+64]);
x0 = Spectrum.periodic(x0);

params.scales  = [1,.5,.25];
t1 = Texture(x0,params);

wc_params.num_levels = 4;
t1.add_constraint(WaveletConstraint(t1,wc_params));
% t1.add_constraint(HistogramConstraint(t1));

mc_params.blocksize = 3;
mc_params.x_dataratio = .25;
mc_params.y_dataratio = .4;
mc_params.num_levels = 4;
mc_params.scales = params.scales;

% t1.add_constraint(MRFConstraint(t1,mc_params));
t1.add_constraint(WaveletAncestryMRFConstraint(t1,mc_params));


t1.run_variational_synthesis()



Analyze.get_tiling(t1.x0,t1.y, 7 , 1);