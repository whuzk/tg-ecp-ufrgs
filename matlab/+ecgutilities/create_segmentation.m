function Result = create_segmentation(Database)
import ecgutilities.*

Records = fieldnames(Database);
for i = 1:numel(Records)
    ECG = Database.(Records{i});
    
    [Fs,Bp,Sc,Leads,Data] = ecgutilities.interpret(ECG);
    for j = 1:Sc
        disp(['processing ' Records{i} '.' Leads{j}]);
        
        % processa o sinal da derivašao
        Filt = ecgfastcode.c_filter_double(Data(:,j),Fs);
        %Rpeaks = ecgfeatures.detect_qrs(Filt,Fs);
        Rpeaks = ecgfeatures.sogari_qrs(Filt,Fs);
        [A,B,R] = ecgutilities.merge_rpeaks(Bp, Rpeaks, Fs);
        
        % salva o resultado
        Result.(Records{i}).(Leads{j}).A = A;
        Result.(Records{i}).(Leads{j}).B = B;
        Result.(Records{i}).(Leads{j}).R = R;
        %{
        figure; hold on;
        plot(Data(:,j));
        stem(R, A*0.2, 'g');
        stem(R, B*0.1, 'r');
        %}
        ecgmath.compute_statistics(A,B)
    end
end