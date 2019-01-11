[file, path] = uigetfile('./*.mph');
abs_path = strcat(path, file);
model = mphopen(abs_path);
% Add export section
data = model.result.export.create('2D_data', 'Data');
data.set('recover', 'pprint'); %use recover setting
data.set('fullprec', 'off');
%use grid location
data.set('location', 'grid');
data.set('gridx2', 'range(0, (R_p*2-0)/511, R_p*2)');
data.set('gridy2', 'range(-R_p*2, (R_p * 4)/511, R_p*2)');

%use
expressions = {'V', 'c_p', 'c_n', 'tds.tflux_c_pz', 'tds.tflux_c_nz'};
desc = {'Potential', 'c_p', 'c_n', 'cpz', 'cnz'};
data.set('expr', expressions);
data.set('descr', desc);

%file
file_root = strsplit(file, '.');
txt_file = strcat(path, strcat(file_root{1}, '_2Ddata.txt'));
data.set('filename', txt_file);

%execute

fprintf('Exporting for file %s ...\n', file);
data.run();