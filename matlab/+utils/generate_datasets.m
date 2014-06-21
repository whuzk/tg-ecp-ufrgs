import utils.*
global leadnames

%basedir = 'C:\physiobank\database\edb\extracted\'; 
basedir = '/Users/diegosogari/Documents/MATLAB/extracted/';
neglect = {'e0403'};

for i = 3:length(leadnames)
    lead = leadnames{i};
    disp(['Processing lead ' lead '...']);
    Rocha.(lead) = get_lead_dataset(basedir, lead, 'Rocha', neglect);
    Mohebbi.(lead) = get_lead_dataset(basedir, lead, 'Mohebbi', neglect);
    Gopalak.(lead) = get_lead_dataset(basedir, lead, 'Gopalak', neglect);
end