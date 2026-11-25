function results = chaudhry_example_8_2()
% CHAUDHRY_EXAMPLE_8_2  Frequency response of a branching pipe system
% by the transfer matrix method.
%
% Reference: M. H. Chaudhry, Applied Hydraulic Transients, 3rd ed.,
% Springer, 2014, Chapter 8, Example 8-2 and Fig. 8-17.
%
% System: reservoir -> pipe 1 -> junction (dead-end branch) -> pipe 2
%         -> oscillating valve discharging to atmosphere.
%
% Method:
%   Field matrix, frictionless pipe (Eq. 8-37)
%       F = [cos(b w)        -j sin(b w)/C ;
%            -j C sin(b w)    cos(b w)     ],  b = L/a,  C = a/(g A)
%   Point matrix, dead-end branch (Eq. 8-87)
%       Pbde = [1  u12/u11 ; 0  1],  u = field matrix of the branch
%   Overall transfer matrix (Fig. 8-17b)
%       U = F2 * Pbde * F1
%   Oscillating valve, single exciter (Eqs. 8-133 to 8-135)
%       q1 = -(2 Ho k / tau0) / (u21 - (2 Ho / Qo) u11)
%       q3 = u11 q1 ,  h3 = u21 q1
%   Normalised response
%       hr = 2|h3|/Ho ,  qr = 2|q3|/Qo
%
% State vector is {q; h}: complex amplitudes of discharge and head
% oscillation. Runs in MATLAB and GNU Octave.
%
% Johann Cruz

    % ---- Data (Example 8-2, Fig. 8-17a) --------------------------------
    g    = 9.80665;   % m/s2, standard gravity
    Ho   = 100;       % m, mean head at the valve
    Qo   = 0.314;     % m3/s, mean discharge
    k    = 0.2;       % amplitude of relative gate opening oscillation
    tau0 = 1.0;       % mean relative gate opening
    Tth  = 3.0;       % s, theoretical period = 4(L1 + L2)/a
    wth  = 2*pi/Tth;  % rad/s, theoretical frequency

    %             L (m)   a (m/s)   D (m)
    pipe1  = struct('L', 500, 'a', 1000, 'D', 1.81);
    pipe2  = struct('L', 250, 'a', 1000, 'D', 1.05);
    branch = struct('L', 200, 'a',  800, 'D', 1.36);

    % ---- Check against the hand calculation at wr = 2.0 ---------------
    r = response(2.0);
    fprintf('\nExample 8-2 at wr = 2.0\n');
    fprintf('%-34s %10s %12s\n', 'Quantity', 'Book', 'This code');
    fprintf('%-34s %10.3f %12.4f\n', 'hr = 2|h3|/Ho',          0.170,  r.hr);
    fprintf('%-34s %10.3f %12.4f\n', 'qr = 2|q3|/Qo',          0.390,  r.qr);
    fprintf('%-34s %10.2f %12.4f\n', 'Phase, head vs gate (deg)',   -102.29, r.phase_h);
    fprintf('%-34s %10.2f %12.4f\n', 'Phase, flow vs gate (deg)',   -12.29,  r.phase_q);
    fprintf('(The book rounds intermediate values in its hand calculation.)\n\n');

    % ---- Frequency sweep ----------------------------------------------
    wr = linspace(0.5, 7.0, 3001);
    hr = zeros(size(wr));
    qr = zeros(size(wr));
    for n = 1:numel(wr)
        s = response(wr(n));
        hr(n) = s.hr;
        qr(n) = s.qr;
    end

    % Resonant peaks of hr
    pk = find(hr(2:end-1) > hr(1:end-2) & hr(2:end-1) > hr(3:end)) + 1;
    fprintf('Resonant peaks of hr:\n');
    fprintf('   wr        w (rad/s)   period (s)   hr\n');
    for n = pk
        w = wr(n)*wth;
        fprintf('  %6.3f    %8.3f    %8.3f    %7.3f\n', wr(n), w, 2*pi/w, hr(n));
    end

    % ---- Plot (compare with Fig. 8-17c) -------------------------------
    figure('Color', 'w');
    plot(wr, hr, 'k-', 'LineWidth', 1.4); hold on;
    plot(wr, qr, 'k--', 'LineWidth', 1.2);
    plot(2.0, r.hr, 'ko', 'MarkerFaceColor', 'w');
    plot(2.0, r.qr, 'ko', 'MarkerFaceColor', 'w');
    grid on; xlim([0.5 7]); ylim([0 0.85]);
    xlabel('\omega_r = \omega / \omega_{th}');
    ylabel('h_r ,  q_r');
    legend('h_r = 2|h_3|/H_o', 'q_r = 2|q_3|/Q_o', 'Location', 'southoutside', ...
           'Orientation', 'horizontal');
    title('Example 8-2: frequency response at the oscillating valve');

    results = struct('wr', wr, 'hr', hr, 'qr', qr, 'check', r);

    % ---- Nested functions ---------------------------------------------
    function out = response(wr_i)
        w   = wr_i*wth;
        F1  = field_matrix(pipe1,  w, g);
        F2  = field_matrix(pipe2,  w, g);
        Fb  = field_matrix(branch, w, g);
        Pbde = [1, Fb(1,2)/Fb(1,1); 0, 1];
        U   = F2*Pbde*F1;
        q1  = -(2*Ho*k/tau0) / (U(2,1) - (2*Ho/Qo)*U(1,1));
        q3  = U(1,1)*q1;
        h3  = U(2,1)*q1;
        out = struct('hr', 2*abs(h3)/Ho, 'qr', 2*abs(q3)/Qo, ...
                     'phase_h', angle(h3)*180/pi, 'phase_q', angle(q3)*180/pi, ...
                     'q1', q1, 'q3', q3, 'h3', h3, 'U', U);
    end
end

function F = field_matrix(p, w, g)
% Frictionless field matrix of a uniform pipe (Chaudhry Eq. 8-37).
    A = pi*p.D^2/4;
    C = p.a/(g*A);          % characteristic impedance, s/m2
    b = p.L/p.a;            % wave travel time, s
    F = [cos(b*w),            -1i*sin(b*w)/C;
         -1i*C*sin(b*w),       cos(b*w)];
end
