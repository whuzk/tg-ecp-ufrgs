function Result = test_networks(Dataset, Network, MethodName)

B = ones(1,10)/10;
F = filter(B,1,Dataset.(MethodName).Features);
D = 2*Dataset.(MethodName).Diagnosis-1;

switch MethodName
    case 'Rocha'
        [T,O] = test_rocha(F, D, Network.RochaNet);
    case 'Mohebbi'
        [T,O] = test_mohebbi(F, D, Network.MohebbiNet);
    case 'Gopalak'
        [T,O] = test_gopalak(F, D, Network.GopalakNets);
    otherwise
        error('invalid method name');
end

Result = ecgmath.compute_statistics(T, O);


function [T,O] = test_rocha(F, D, Network)
A = sim(Network, F')';
O = A(:,1) > 0 | A(:,2) > 0;
T = D(:,1) > 0 | D(:,2) > 0;

function [T,O] = test_mohebbi(F, D, Network)
A = sim(Network, F')';
O = A(:,1) > 0 & A(:,2) < 0;
T = D(:,1) > 0 & D(:,2) < 0;

function [T,O] = test_gopalak(F, D, Networks)
k = length(Networks);
m = size(F,1);
decision = zeros(m,1);
for i = 1:k
    A = sim(Networks{i}, F')';
    C = A(:,1) > 0 | A(:,2) > 0;
    decision = decision + double(C);
end
O = decision > k/2;
T = D(:,1) > 0 | D(:,2) > 0;
