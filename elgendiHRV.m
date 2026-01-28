%% Inicijalizacija

clc; clearvars;
set(0, 'DefaultFigureVisible', 'off');

% ucitavanje .fig filea
g = load('filepath.fig','-mat'); % dodatak vlastitog patha datoteke
x2 = g.hgS_070000.children(3,1).children(1,1).properties.XData;
y2 = g.hgS_070000.children(3,1).children(1,1).properties.YData;

Fs = 364; 
t = x2;
shift = 0.1;

%% Filtriranje
bpfilt = designfilt('bandpassiir','FilterOrder',2, 'HalfPowerFrequency1',0.5,'HalfPowerFrequency2',8, 'SampleRate',Fs);

y_filt = filtfilt(bpfilt, y2); 
y = (y_filt - min(y_filt)) ./ (max(y_filt) - min(y_filt)); % normalizacija

%% Generiranje blokova
w1 = round(0.11 * Fs); 
w2 = round(0.65 * Fs); 

map = movmean(y, w1);
mab = movmean(y, w2);

beta = 0.001;
alpha = beta * mean(y);
blockofinterest = map > (mab + alpha);

%% Elgendi peak Detection (anti double count)
i = 1; k = 1; locs = [];
refractory_samples = round(0.40 * Fs); % 400ms zabrane detekcije

while i < length(blockofinterest)
    if blockofinterest(i) == 1
        start_idx = i;
        while i < length(blockofinterest) && blockofinterest(i) == 1
            i = i + 1;
        end
        end_idx = i;
        
        [~, pki] = max(y_filt(start_idx:end_idx));
        current_loc = (start_idx + pki - 1) / Fs;
        
        % provjera je li novi vrh predaleko od prošlog
        if k == 1 || (current_loc - locs(k-1)) > 0.40 % 400 milisekundi
            locs(k) = current_loc;
            k = k + 1;
            i = i + refractory_samples; % Preskoči refraktorni period
        end
    else
        i = i + 1;
    end
end

%% HRV filtriranje (augmented)

HR_ibi = diff(locs);

%postavljanje prihvatljivog opsega za HR
valid_ibi = HR_ibi(HR_ibi > 0.45 & HR_ibi < 1.3);

% medijan filter
med_ibi = median(valid_ibi);
HR_ibi_final = valid_ibi(valid_ibi > med_ibi*0.85 & valid_ibi < med_ibi*1.15);

% filter razlika
rr_ms = HR_ibi_final * 1000;
diff_rr = diff(rr_ms);
valid_diffs = abs(diff_rr) < 80; 

rr_ms_clean = rr_ms([true, valid_diffs]); 
diff_rr_clean = diff(rr_ms_clean);

%% Ispis proračunatih HRV varijanta u vremenskoj domeni
SDNN = std(rr_ms_clean);
RMSSD = sqrt(mean(diff_rr_clean.^2));
pNN50 = (sum(abs(diff_rr_clean) > 50) / length(diff_rr_clean)) * 100;

fprintf('\nREZULTATI:\n');
fprintf('Broj otkucaja: %d \n', length(rr_ms_clean));
fprintf('Prosječni HR: %.1f bpm\n', 60000 / mean(rr_ms_clean));
fprintf('SDNN: %.2f ms \n', SDNN);
fprintf('RMSSD: %.2f ms \n', RMSSD);
fprintf('pNN50: %.2f %% \n', pNN50);

%kraj