function convert_all_databases
% converte as bases de dados de ECGs

% converte a base European completa
fprintf('Converting ST-T complete database...\n');
desc_dir = 'C:\ecg\edb\hea\';
data_dir = 'C:\ecg\edb\sgnl\';
atr_dir = 'C:\ecg\edb\atr\';
save_file = 'C:\ecg\edb.mat';
utilities.convert_database(...
    desc_dir, data_dir, atr_dir, '', '', '', '', save_file);
pause(1);

% converte a base LT-ST completa
fprintf('Converting LT-ST complete database...\n');
desc_dir = 'C:\ecg\ltstdb\hea\';
data_dir = 'C:\ecg\ltstdb\sgnl\';
atr_dir = 'C:\ecg\ltstdb\atr\';
sta_dir = 'C:\ecg\ltstdb\sta\';
stb_dir = 'C:\ecg\ltstdb\stb\';
stc_dir = 'C:\ecg\ltstdb\stc\';
save_file = 'C:\ecg\ltstdb2.mat';
utilities.convert_database(...
    desc_dir, data_dir, atr_dir, '', sta_dir, stb_dir, stc_dir, save_file);
pause(1);

% converte a base QT completa
fprintf('Converting QT complete database...\n');
desc_dir = 'C:\ecg\qtdb\hea\';
data_dir = 'C:\ecg\qtdb\sgnl\';
atr_dir = 'C:\ecg\qtdb\atr\';
pu_dir = 'C:\ecg\qtdb\pu\';
save_file = 'C:\ecg\qtdb.mat';
utilities.convert_database(...
    desc_dir, data_dir, atr_dir, pu_dir, '', '', '', save_file);
pause(1);

% converte a base MIT-BIH completa
fprintf('Converting MIT-BIH complete database...\n');
db_dir = 'C:\ecg\mitdb\';
save_file = 'C:\ecg\mitdb.mat';
utilities.convert_physiobank(db_dir, save_file);
