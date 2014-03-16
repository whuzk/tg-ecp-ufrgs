function [R,Th] = tompkins_adapted(Signal, Fs)

%initializations
Sthresh = 0;
Nthresh = 0;
Dthresh = 0;
achtung = false;
LASTPEAK = 0;
N = length(Signal);
M = fix(N/Fs*5);
R = zeros(M,1);
Th = zeros(N,1);
FOsize = 20;

% algorithm
k = 0;
i = 2;
while k <= M && i <= N-FOsize
    % get signal amplitude at current position
    P = Signal(i);
    
    % update global threshold according to the signal and noise thresholds
    Gthresh = Nthresh + 0.25*(Sthresh-Nthresh);
    
    % change state of the algorithm depending on whether the signal is
    % above the treshold or not
    if P > Gthresh
        if ~achtung
            % algorithm is entering the state of "active search" (when a
            % new QRS peak must be detected or refused)
            achtung = true;
            %i, pause;
        end
    elseif achtung
        % algorithm is leaving the state of active search
        achtung = false;
        % reset the last QRS peak (we know that the signal is in the range
        % of positive real values, so we can use zero as a special marker)
        LASTPEAK = 0;
        %i, pause;
    end
    
    % check if the current point is a positive peak or the verge of one
    a = Signal(i-1);
    c = Signal(i+1);
    if (P-a) > 0 && (c-P) <= 0 && abs(c-2*P+a) > 1E-9
        if ~achtung
            % peak has been considered as noise peak
            Nthresh = 0.125*P + 0.875*Nthresh;
            % give a decrease on the signal threshold, to quickly recover
            % from possible (and unusual) high peaks in the signal
            Sthresh = Sthresh - 0.025*(Sthresh-Nthresh);
            %i, pause
        elseif P-Gthresh < Dthresh
            % peak has been considered as artifact peak
            Dthresh = 0.5*Dthresh + 0.1*(P-Gthresh);
            %i, pause
        elseif P > Signal(i+FOsize)
            % peak has been considered as signal peak
            Sthresh = 0.125*P + 0.875*Sthresh;
            if P > LASTPEAK
                % peak is a new high for the current QRS
                Dthresh = 0.2*(P-Gthresh);
                if LASTPEAK == 0
                    % peak is the first for the current QRS
                    k = k + 1;
                end
                LASTPEAK = P;
                R(k) = i;
                %i, pause
            end
        end
    end
    Th(i) = Gthresh;
    i = i + 1;
end
R(k+1:end) = [];