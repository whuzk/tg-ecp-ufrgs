function convert_all_databases
% converte as bases de dados de ECGs
import ecgutilities.*
%{
% converte a base European completa
fprintf('Converting ST-T complete database...\n');
db_dir = 'I:\AppData\physiobank\database\edb';
save_file = 'C:\ecg\edb.mat';
convert_physiobank_db(db_dir, save_file);
pause(1);

% converte a base LT-ST completa
fprintf('Converting LT-ST complete database...\n');
db_dir = 'I:\AppData\physiobank\database\ltstdb';
save_file = 'C:\ecg\ltstdb.mat';
convert_physiobank_db(db_dir, save_file);
pause(1);
%}
% converte a base QT completa
fprintf('Converting QT complete database...\n');
db_dir = 'I:\AppData\physiobank\database\qtdb';
save_file = 'C:\ecg\qtdb.mat';
convert_physiobank_db(db_dir, save_file);
pause(1);
%{
% converte a base MIT-BIH completa
fprintf('Converting MIT-BIH complete database...\n');
db_dir = 'I:\AppData\physiobank\database\mitdb\';
save_file = 'C:\ecg\mitdb.mat';
convert_physiobank_db(db_dir, save_file);
%}