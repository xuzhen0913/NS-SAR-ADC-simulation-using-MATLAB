clc;
clear;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% 定义SAR ADC的基本参数 %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

OSR = 8;              % 过采样率
N   = 12;             % ADC resolution (bit)
fs  = 1e6/OSR;        % 采样频率 (Hz)
Vref = 1;             % 参考电压 (V)
Vos  = 0.5e-3;        % comparator offset (V)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% 定义环境参数 %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

k = 1.38e-23;
T = 300;
sigma_buffer = 200e-6;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% 定义CDAC的基本参数 %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Cu = (1e-12)/(2^N);
C_arr_p = [2.^[N-1:-1:0],1] .* Cu;
C_arr_n = [2.^[N-1:-1:0],1] .* Cu;
C_tot_p = sum(C_arr_p);
C_tot_n = sum(C_arr_n);

sigma_v_p = sqrt(k*T / C_tot_p);
sigma_v_n = sqrt(k*T / C_tot_n);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% 定义输入信号参数 %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

num = 2^10;
Vfs = Vref;
fin = (13/num)*fs;
ts  = 1/fs;
t   = (0:num-1)' * ts;

Vin_p = Vref/2 + (Vfs/2)*sin(2*pi*fin*t + pi/3);
Vin_n = Vref/2 - (Vfs/2)*sin(2*pi*fin*t + pi/3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% <<< NS-SAR ADD :  EF >>> %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% -------- EF : FIR ----------
% H_EF(z) = (1 - z^-1)
b_ef = [1 -1];
a_ef = 1;

zi_ef  = zeros(length(b_ef)-1,1);
e_prev = 0;   % quantization error memory

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% ADC工作过程 %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Dout_10 = zeros(num,1);

for i = 1:num

    code_osr = zeros(OSR,1);

    for m = 1:OSR

        % ===== 噪声 =====
        vn_p = sigma_v_p * randn;
        vn_n = sigma_v_n * randn;
        vn_buffer_p = sigma_buffer * randn;
        vn_buffer_n = sigma_buffer * randn;

        Vin_p_osr = Vin_p(i) + vn_p + vn_buffer_p;
        Vin_n_osr = Vin_n(i) + vn_n + vn_buffer_n;

        Vin_diff_osr = Vin_p_osr - Vin_n_osr;

        % ===== <<< EF : error path >>> =====
        [x_ef, zi_ef] = filter(b_ef, a_ef, e_prev, zi_ef);

        % ===== 量化器等效输入 =====
        Vin_eff = Vin_diff_osr-x_ef;

        % ===== SAR conversion =====
        V_p_bottom = zeros(1,N+1);
        V_n_bottom = Vref*ones(1,N+1);

        Dout_tmp = zeros(1,N);

        for j = 1:N
            V_p_bottom(j) = Vref;
            V_n_bottom(j) = 0;

            Vtop_p = (C_arr_p * V_p_bottom') / C_tot_p - Vin_eff/2;
            Vtop_n = (C_arr_n * V_n_bottom') / C_tot_n + Vin_eff/2;

            if (Vtop_p - Vtop_n + Vos) < 0
                Dout_tmp(j) = 1;
            else
                Dout_tmp(j) = 0;
                V_p_bottom(j) = 0;
                V_n_bottom(j) = Vref;
            end
        end

        code = Dout_tmp * 2.^[(N-1):-1:0]';
        code_osr(m) = code;

        % ===== 更新量化误差 =====
        Vq = (code/2^N - 0.5) * 2*Vref;
        e_prev = Vq-Vin_eff;

    end

    Dout_10(i) = mean(code_osr) / 2^N;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% 后处理 %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Vin_diff  = Vin_p - Vin_n;
Vout_diff = (Dout_10 - 0.5) * (2*Vref);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% 作图 %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;
plot(t, Vin_diff, 'b'); hold on;
stairs(t, Vout_diff, 'r');
grid on;
xlabel('Time (s)');
ylabel('Voltage (V)');
legend('Input','Output');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% 频域性能 %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x = Vout_diff;

X = fft(x);

%X = X(1:num/2);
%P = abs(X).^2;
%f = (0:num/2-1)'*fs/num;

P2 = abs(X/num);
P1 = P2(1:num/2+1);
Ppow = P1.^2;
P1(2:end-1) = 2*P1(2:end-1);
f = fs/num*(0:(num/2));

k_sig = 14;
P_signal = sum(P1(k_sig-1:k_sig+1).^2);
P_noise  = sum(P1(2:end).^2) - P_signal;

SNR  = 10*log10(P_signal/P_noise);
ENOB = (SNR - 1.76)/6.02;

fprintf('SNR  = %.2f dB\n', SNR);
fprintf('ENOB = %.2f bits\n', ENOB);

figure;
plot(f/1e3,10*log10(Ppow / max(Ppow)));
xlabel('Frequency (kHz)');
grid on;
xlabel('Frequency (kHz)');
ylabel('Magnitude (dBFS)');