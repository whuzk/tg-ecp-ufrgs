function [QRS,QRSreal,RR] = tompkins_original(Signal, Fs)

%lahead = round(1:2*Fs);
%noiseUpdate = .125;
%signalUpdate = .125;
%signalSBupdate = .25;
%rrHighPerc = 1.16;
%rrLowPerc = .92;
%rrMissLim = 1.66;
peakDur = ceil(75/1000*Fs);
QRSupdate = .475;
Resolution = 12;
lenRRavg = 8;

%RRmiss = 0;
%searchback = 0;
%waiting = 0;
%lastQRS = 1;
%waitForNextPeak = 0;
%npk = min(Signal(lahead));
%spk = max(Signal(lahead));
%blankPer = round(.2*Fs);
%maxder = 0;
%range = (-peakDur+1):0;

thresh1 = 2/2^12;
thresh2 = thresh1/2;
qpkcnt = 0;
cnt = 0;
sbpeak = 0;
initBlank = 0;
initMax = 0;
preBlankCnt = 0;
sbCount = round(1.5*Fs);
tempPeak = 0;
WINDOW_WIDTH = peakDur;
PRE_BLANK = round(2*peakDur);
pos = 1;
v2 = 0;
rsetBufPos = 1;
noisePos = 1;
timeSinceMax = 0;

RR = zeros(size(Signal));
RR(1:lenRRavg) = Fs;
rsetBuf = zeros(1,lenRRavg);
noisePeakBuf = zeros(1,lenRRavg);
QRS = zeros(size(Signal));
QRSreal= zeros(size(Signal));
QRSpeakBuf = zeros(size(Signal));
thresh1ar = zeros(size(Signal));
thresh1ar(1:pos) = thresh1;
thresh2ar = zeros(size(Signal));
thresh2ar(1:pos) = thresh2;

dMax = 0;
pkTmp = 1;
pkTmp2 = 1;
pkPos = 1;
pkCnt = 1;
qmedian = 0;
nmedian = 0;

while pos <= length(Signal)

    v = Signal(pos);

    if (timeSinceMax > 0) 
        timeSinceMax  =timeSinceMax+1;
    end
    
    pk = 0;
    
    if v > v2 && v > dMax
        dMax = v;
        pkTmp2 = pos;
        if dMax > 2/2^Resolution 
            timeSinceMax = 1; 
        end
    elseif v < (dMax/2)
        pk = dMax;                
        dMax = 0;
        timeSinceMax = 0;
    elseif timeSinceMax > round(.095*Fs)
        pk = dMax;
        dMax = 0;
        timeSinceMax = 0;
    end

    newPeak = 0;

    if (pk ~= 0 && preBlankCnt==0)
        tempPeak = pk;
        pkTmp = pkTmp2;
        preBlankCnt = PRE_BLANK;
    elseif pk==0 && preBlankCnt ~= 0
        preBlankCnt = preBlankCnt-1;
        if preBlankCnt == 0 
            newPeak = tempPeak; 
            pkPos = pkTmp;
        end
    elseif pk ~= 0
        if pk > tempPeak
            tempPeak = pk;
            pkTmp = pkTmp2;
            preBlankCnt = PRE_BLANK;
        else
            preBlankCnt = preBlankCnt-1;
            if preBlankCnt == 0
                newPeak = tempPeak;
                pkTmp = pkTmp2;
            end
        end
    end
    
    if qpkcnt < 9
        cnt = cnt+1;
        if newPeak > 0 
            cnt = WINDOW_WIDTH;
        end
        initBlank = initBlank+1;
        if initBlank == Fs || newPeak > 0
            if newPeak > thresh1
                initBlank = 0;
                QRSpeakBuf(pkCnt) = newPeak;
                QRS(pkCnt) = pos;
                QRSreal(pkCnt) = pkPos;
                pkCnt = pkCnt+1;
                initMax = 0;
                qpkcnt = qpkcnt+1;
                %if qpkcnt==lenRRavg
                qmedian = median(QRSpeakBuf(1:qpkcnt));
                nmedian = 0;
                %rrmedian = Fs;
                %sbcount = round(Fs*(1.500+.150));
                %dmed = qmedian-nmedian;
                thresh1 = calcThresh(qmedian, nmedian, QRSupdate);
                thresh2 = thresh1/2;
            else
                noisePeakBuf(noisePos) = newPeak;
                noisePos = noisePos+1;
                if noisePos > 8
                    noisePos=1;
                end                        
                nmedian = median(noisePeakBuf);
                thresh1 = calcThresh(qmedian, nmedian, QRSupdate);
                if qmedian == 0
                    thresh1 = mean(noisePeakBuf)*(1+QRSupdate);
                end
                thresh2 = thresh1/2;
            end
            %end
        end
        if newPeak > initMax
            initMax = newPeak;
            pkPos = pkTmp;
        end
    else
        cnt = cnt+1;
        if newPeak > 0
            if newPeak > thresh1
                %qrsDelay = signal.totalDelay;
                QRSpeakBuf(pkCnt) = newPeak;
                QRS(pkCnt) = pos;%-qrsDelay;
                QRSreal(pkCnt) = pkTmp;
                qmedian = median(QRSpeakBuf(pkCnt-7:pkCnt));
                thresh1 = calcThresh(qmedian, nmedian, QRSupdate);
                RR(pkCnt) = cnt;
                pkCnt = pkCnt+1;
                %rrmedian = median(RR(pkCnt-7:pkCnt));
                %sbcount = rrmedian + rrmedian/2 + peakDur;
                cnt = peakDur;
                sbpeak = 0;
                %lastmax = maxder;
                %maxder = 0;
                initBlank = 0;
                initMax = 0;
                rsetBufPos = 1;
            else
                noisePeakBuf(noisePos) = newPeak;
                noisePos = noisePos+1;
                if noisePos > lenRRavg
                    noisePos = 1;
                end                        
                nmedian = median(noisePeakBuf);
                thresh1 = calcThresh(qmedian, nmedian, QRSupdate);
                thresh2 = thresh1/2;

                if (newPeak > sbpeak) && (cnt - peakDur >= round(Fs*.36))
                    sbpeak = newPeak;
                    sbloc = cnt - peakDur;
                end
            end
        end
        if cnt > sbCount && sbpeak > thresh2
            qrsDelay = cnt - sbloc;
            cnt = cnt-sbloc;
            %qrsDelay = qrsDelay + signal.totalDelay;

            QRSpeakBuf(pkCnt) = sbpeak;
            QRS(pkCnt) = pos-qrsDelay;
            QRSreal(pkCnt) = pkTmp;

            qmedian = median(QRSpeakBuf(pkCnt-lenRRavg:pkCnt-1));
            thresh1 = calcThresh(qmedian, nmedian, QRSupdate);
            thresh2 = thresh1/2;
            RR(pkCnt) = sbloc;
            pkCnt = pkCnt+1;
            %rrmedian = median(RR(pkCnt-lenRRavg:pkCnt-1));
            %sbcount = rrmedian + rrmedian/2 +peakDur;                   
            
            sbpeak = 0;
            %lastmax = maxder;
            %maxder = 0;
            initBlank = 0;
            initMax = 0;
            rsetBufPos = 1;
        end
    end
    if qpkcnt == lenRRavg
        initBlank = initBlank+1;
        if initBlank == Fs
            initBlank = 0;
            rsetBuf(rsetBufPos) = initMax;
            rsetBufPos = rsetBufPos+1;
            if rsetBufPos > lenRRavg
                rsetBufPos=1;
            end
            initMax = 0;
            if rsetBufPos == lenRRavg
                for q=1:lenRRavg
                    QRSpeakBuf(q) = rsetBuf(q);
                    noisePeakBuf(q) = 0;
                end
                qmedian = median(rsetBuf);
                nmedian = 0;
                %rrmedian = Fs;
                %sbcount = round(Fs*(1.5+.15));
                thresh1 = calcThresh(qmedian, nmedian, QRSupdate);
                thresh2 = thresh1/2;
                initBlank = 0;
                initMax = 0;
                %rsetCount = 0;
                rsetBufPos = 1;
                sbpeak = 0;
            end
        end
        if newPeak > initMax
            initMax = newPeak;
        end
    end
%             [pktmp,idx] = max(Signal(pos+range));
%             if pktmp ~= pk && abs(idx - length(range)/2) <3
%             %if pktmp > Signal(pos+range(1)) && pktmp > Signal(pos+range(end)) && pktmp ~= pk
%                 pk = pktmp;
%                 if ~searchback
%                     if pk >= thresh1
%                         spk = (signalUpdate)*pk + (1-signalUpdate)*spk;
%                     else
%                         npk = noiseUpdate*pk +(1-noiseUpdate)*npk;
%                     end
%                 else
%                     if pk >= thresh2
%                         spk = (signalSBupdate)*pk + (1-signalSBupdate)*spk;
%                     else
%                         npk = noiseUpdate*pk +(1-noiseUpdate)*npk;
%                     end
%                 end
%                 thresh1 = npk + (spk-npk)/4;
%                 thresh2 = thresh1/2;
%             end
%             if ~searchback threshtmp = thresh1;
%             else threshtmp = thresh2;
%             end
%             if v <= spk/2 && v2 > spk/2 && pk > threshtmp && waiting == 0
% %                 if ~searchback
% %                     spk = (signalUpdate)*pk + (1-signalUpdate)*spk;
% %                 else
% %                     spk = (signalSBpdate)*pk + (1-signalSBUpdate)*spk;
% %                 end
%                 QRS = [QRS, pos];
%                 lastQRS = pos;
%                 waitForNextPeak = 0;
%                 if length(QRS) > 1
%                     RR = [RR, QRS(end)-QRS(end-1)];
%                     rRange = min(length(RR), lenRRavg)-1;
%                     RRavg1 = mean(RR(end-rRange:end));
%                     RRavg2 = RRavg1;
%                     if rRange >= 7
%                         % start updating searchback
%                         tmpPos = intersect(find(RR <= RRhigh), find(RR >= RRlow));
%                         tmpRange = min(length(tmpPos),lenRRavg)-1;
%                         RRavg2 = mean(RR(end-tmpRange:end));
%                         RRhigh = rrHighPerc*RRavg2;
%                         RRlow = rrLowPerc*RRavg2;
%                         RRmiss = rrMissLim*RRavg2;
%                     else
%                         RRhigh = rrHighPerc*RRavg2;
%                         RRlow = rrLowPerc*RRavg2;
%                         RRmiss = rrMissLim*RRavg2;
%                     end
% %                     RRavg1(RRavg1Pos) = RR(end) - RR(end-1);
% %                     RRavg1Pos = RRavg1Pos+1;
% %                     if RRavg1Pos >lenRRavg RRavg1Pos = 1; end
%                 end
%                 
%                 waiting = blankPer;
%                 thresh1 = npk + (spk-npk)/4;
%                 thresh2 = thresh1/2;
%                 searchback = 0;
%             end
%             if waiting > 0 
%                 waiting = waiting-1; end
%             
%             % check if searchback is necessary
%             if RRmiss ~= 0
%                 if (pos-lastQRS) >= RRmiss
%                     if searchback == 0 && waitForNextPeak == 0
%                         pos2 = pos;
%                         pos = lastQRS;
%                         v2 = Signal(pos-1);
%                         searchback = 1;
%                         waiting = blankPer;
%                         continue;
%                     else
%                         waitForNextPeak = 1;
%                     end
%                 end
%             end
%             if searchback && pos == pos2
%                 searchback = 0;
%             end
%             
%             v2 = Signal(pos);
    thresh1ar(pos) = thresh1;
    thresh2ar(pos) = thresh2;
    v2 = v;
    pos = pos+1;
end
QRS(pkCnt:end) = [];
QRSreal(pkCnt:end) = [];
RR(pkCnt:end) = [];

function thrsh = calcThresh(qmedian, nmedian, TH)
thrsh = nmedian + TH*(qmedian - nmedian);