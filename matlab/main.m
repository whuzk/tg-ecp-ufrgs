% main.m
%   Programa principal. Nota: execute o arquivo init.m primeiro.
%
import ecgutilities.*
import ecgtest.*
%close all;

%% Criaçao do arquivo de teste do algoritmo de segmentaçao
Segment = create_segmentation(MITDB);
%fprintf('Saving data to file...\n');
%save_filepath = './resources/segment_V3.mat';
%save(save_filepath, '-struct', 'Segment');
%fprintf('Segments saved to %s\n', save_filepath);
%Segment = load('./resources/segment.mat');
%test_segmentation(Segment)

%% Criaçao do arquivo de batidas
%Beatset = create_beatset(EDB);
%fprintf('Saving data to file...\n');
%save_filepath = './resources/beatset.mat';
%save(save_filepath, '-struct', 'Beatset');
%fprintf('Beatset saved to %s\n', save_filepath);
%Beatset = load('./resources/beatset.mat');

%% Criaçao do arquivo de caracteristicas
%Dataset = create_dataset(Beatset, LeadName);
%fprintf('Saving data to file...\n');
%save_filepath = './resources/dataset.mat';
%save(save_filepath, '-struct', 'Dataset');
%fprintf('Dataset saved to %s\n', save_filepath);
%Dataset = load('./resources/dataset.mat');

%% Criaçao da rede neural do metodo de Rocha
%RochaNet = train_network(Dataset.Rocha);
%fprintf('Saving data to file...\n');
%save_filepath = './resources/networks.mat';
%save(save_filepath, '-struct', 'RochaNet', '-append');
%fprintf('Network saved to %s\n', save_filepath);

%% Criaçao da rede neural do metodo de Mohebbi
%MohebbiNet = train_network(Dataset.Mohebbi);
%fprintf('Saving data to file...\n');
%save_filepath = './resources/networks.mat';
%save(save_filepath, '-struct', 'MohebbiNet', '-append');
%fprintf('Network saved to %s\n', save_filepath);

%% Criaçao das redes neurais do metodo de Gopalakrishnan
%for i = 1:5
%   GopalakNets{1} = train_network(Dataset.Gopalak);
%end
%fprintf('Saving data to file...\n');
%save_filepath = './resources/networks.mat';
%save(save_filepath, '-struct', 'GopalakNets', '-append');
%fprintf('Networks saved to %s\n', save_filepath);

%% Teste das redes neurais
%Networks = load('../resources/networks.mat');
%Stats(:,1) = test_networks(Dataset, Networks, 'Rocha');
%Stats(:,2) = test_networks(Dataset, Networks, 'Mohebbi');
%Stats(:,3) = test_networks(Dataset, Networks, 'Gopalak');
%Result = dataset([{Stats} Methods], 'ObsNames', Measures);
