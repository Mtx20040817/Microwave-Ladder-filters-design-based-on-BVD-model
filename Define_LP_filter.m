
%   Function used for the generation of low pass chebyshev filtering
%   function
function [S11_Zeros,E_Zeros,LP_w,e1,er] = Define_LP_filter (fL,fH,Tz,N,RL,omega)
%   Input
%   fL: lower band
%   fH: Upper band
%   Tz: Predefined finite transmission zeros, ascending order
%   RL: Return loss in dB
%   w: independent variable w,check this
%   N: degree
%   omega: range of the horizontal axis

%   Output
%   S11 zeros
%   S12 zeros is the predefined Tzs
%   E zeros

%   Input validation
    % if fL >= fH, error('fL must be < fH'); end
    % if any(diff(Tz) <= 0), error('Tz must be strictly ascending'); end
    % if length(Tz) > N, error('Number of finite TZs cannot exceed N'); end
    syms w  %real value

    f0 = sqrt(fL*fH);
    FBW = (fH-fL)/f0;
    %   LP_w: Tzs after the bandpass-lowpass transformation
    LP_w = zeros(1, length(Tz));  
    %   Bandpass transformation    
    % for i = 1:length(Tz)
    %     LP_w(i) = (Tz(i)/f0-f0/Tz(i))/FBW;
    % end
    
    %test
    for i = 1:length(Tz)
         LP_w(i) = Tz(i);
    end
    
    w0 = sqrt(w^2-1);
    
    %   Recursive generation of function
    u = 1;  %initial u
    v = 0;  %initial v

    %   finite TZ
    for j = 1:length(LP_w)
        u_new = w*u - u/LP_w(j) + w0*v*sqrt(1-1/(LP_w(j)^2));
        v_new = w*v - v/LP_w(j) + w0*u*sqrt(1-1/(LP_w(j)^2));
        u = u_new;
        v = v_new;

    end
    %   Tz at INF
    if (N-length(LP_w)) > 0
        for k = 1:(N-length(LP_w))
            u_new = w*u+w0*v;
            v_new = w*v+w0*u;
            u = u_new;
            v = v_new;
        end 
    end
    uc = sym2poly(u);
    S11_Zeros = roots(uc);  %roots are real value

    %   Generate the function of P and F based on Given RL
    %   Calculate the equalripple factor, check this part
    x = 1i;
    P = 1;  % initial p
    LP_w = 1i*LP_w;
    S11_Zeros = 1i*S11_Zeros;
    %   Generate the polynomial P 
    for i = 1:length(LP_w)
        P = P*(x-LP_w(i));
    end
    
    %   Generate polynomial F
    F = 1;  %initial F
    for j = 1:length(S11_Zeros)
        F = F*(x-S11_Zeros(j));
    end
    e = abs(P/F)*(1/sqrt(10^(RL/10)-1));
    %e=P/(F*sqrt(10^(RL/10)-1)); %ripple factor: ratio of the e and er
    e1_square = 1 + (1/(10^(RL/10)-1))*(abs(P/F))^2;
    e1 = sqrt(e1_square);
    %e1 = abs(sqrt(1+e^2));
    er  = e1/(sqrt(e1^2 - 1));
    %er = abs(sqrt(1+e^2)/e);

    %   construct the polynomial based on complex s
    syms s
    F = 1;  %initial F
    for j = 1:length(S11_Zeros)
        F = F*(s-S11_Zeros(j));
    end
    P = 1;  % initial p
    for i = 1:length(LP_w)
        P = P*(s-LP_w(i));
    end
    
    %Fc = sym2poly(F);
    %Pc = sym2poly(P);
    Internal = P/abs(e1) - 1i*F/abs(er);
    Ic = sym2poly(Internal);
    E_Zeros = roots(Ic);
    for k = 1:length(E_Zeros)
        if real(E_Zeros(k)) > 0
            E_Zeros(k) = -conj(E_Zeros(k));
        end
    end

    % The module to plot the polynomial
    count = 1;
    %   omega is the input range based on real value
    for x_range = omega
        x=x_range*1i;    %complex variable

        P = 1;  % initial p
        for i = 1:length(LP_w)
            P = P*(x-LP_w(i));
        end
        F = 1;  %initial F
        for j = 1:length(S11_Zeros)
            F = F*(x-S11_Zeros(j));
        end
        % E = 1;  %initial E
        % for k = 1:length(I_Zeros)
        %     E = E*(x-I_Zeros(k));

        %end
        C=F/P;
        buf=1/(1+(e^2)*(C^2));  %ignore the k in here?
        S21(count)=20*log10((sqrt(buf)));
        S11(count)=20*log10((sqrt(1-buf)));
        count=count+1;

    end
    
    xx= omega;
    plot(xx,S21,'r','linewidth',2);hold on;
    plot(xx,S11,'b','linewidth',2);
    %hold off;
    grid on;
    set(gca,'linewidth',2)
    xlabel('LOWPASS PROTOTYPE FREQUENCY (rad/dec)','fontsize',14)
    ylabel('RETURN LOSS (dB)','fontsize',14)
    legend('S21','S11');




end






