[file, path] = uigetfile('./*.mph');
abs_path = strcat(path, file);
model = mphopen(abs_path);
% Add export section
data = model.result.export.create('1D_data', 'Data');
data.set('data', 'cln2'); %use line cln2
data.set('recover', 'pprint'); %use recover setting
data.set('fullprec', 'off');
%use grid location
data.set('location', 'regulargrid');
data.set('regulargridx1', 1024);

%use
expressions = {'V', 'c_p', 'c_n', 'tds.tflux_c_pz', 'tds.tflux_c_nz'};
desc = {'Potential', 'c_p', 'c_n', 'cpz', 'cnz'};
data.set('expr', expressions);
data.set('descr', desc);

%file
file_root = strsplit(file, '.');
txt_file = strcat(path, strcat(file_root{1}, '_1Ddata.txt'));
data.set('filename', txt_file);

%execute

fprintf('Exporting for file %s ...\n', file);
data.run();

% Export sigma and export file
int_surf = model.result.numerical.create('int2', 'IntSurface');
int_surf.set('intvolume', 'on');
int_surf.selection.set([1, 2]);
int_surf.setIndex('expr', '(z_n * c_n + z_p * c_p) * N_A_const * e_const / (pi * L ^2)', 0);
res = int_surf.computeResult(); cc = res(1); sigma = cc{1};

% Vg = [0.001, 0.025, 0.05, 0.1, 0.15, 0.2, 0.25, 0.35, 0.45, 0.55, 0.65];
Vg = [0.001, 0.025 : 0.025: 0.30];
sigma_file = strcat(path, 'sigma_', file_root{1}, '.txt');
dlmwrite(sigma_file, [Vg', sigma], 'delimiter',' ','precision',5)
