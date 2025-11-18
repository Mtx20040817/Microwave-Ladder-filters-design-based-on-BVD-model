function [B,Jr,Bn1,Bn2] = Parameter_extraction(Yin,Ji_square,Tz,GL)
%   EXTRACT_LADDER_LOWPASS_PROTOTYPE
%   Implements the nodal parameter extraction for symmetric AW ladder filters
%   as described in Giménez et al., IEEE Access 2018
%
%   Inputs:
%   tz_norm   : Vector of normalized transmission zeros [omega1, ..., omegaN]
%               (positive = series, negative = shunt)
%   Yin_num_c,Yin_denum_c:  input impedance
%   Yin:    syms Yin, type is the symbolic expression
%   Ji_square:  input design parameter
%   RL_dB     : Return loss in dB (e.g., 22)
%   GL: prescribed load conductance
%   Deal with the case when the number is too large

%   Outputs:
%   J         : Admittance inverters [J1, Jr1, J2, Jr2, ..., Jn+1]
%   B         : NRN susceptances (FIRs) [B1, B2, ..., Bn]
%   b         : Resonator frequency offsets [b1, ..., bn] = -tz_norm
%   theta11   : Input phase shift
%   yin_poly  : Final input admittance polynomial (for verification)
    syms s
    Ji2 = Ji_square;    %square value

    B = [];
    Jr = [];    %square value
    count = 0;
    %   extraction of dangling resonator
    %   calculate the input phase first


    for i = 1:(length(Tz)-1)

        %follow the procedure from amari's paper
        [Yin_num,Yin_denum] = numden(Yin);  %   sym polynomials, with variable s
        %   bug with large number
        Y_inverse = Yin_denum/Yin_num;
        Yin_nc = double(sym2poly(Yin_num)); 
        Yin_dc = double(sym2poly(Yin_denum));   
        [res, poles,k] = residue(Yin_dc,Yin_nc);  % for Y_inverse
        % disp(res)
        % disp(k)
        % disp(poles)
        %   residue and pole's index vertification, debug this and
        %   optimizee it
        index = 1;
        for m = 1:length(poles)
            if (abs(real(poles(m))) - 0) < 1e-6
                index = m;
            end
        end
        % disp(index)
        Jri2 = Ji2*res(index);  % the inverter for resonator branch
        Jri2 = real(Jri2);
        %disp(Jri2)  
        Jr = [Jr Jri2];
        % modify based on location
        res(index) = [];
        %disp(res)
        poles(index) = [];
        %disp(poles)
        %poly1 = s - poles(2*i-1);  %check the variable, debug this part
        % check this line
        [Yin_nc,Yin_dc] = residue(res,poles,k);
        Yin_num = poly2sym(Yin_nc,s);
        Yin_denum = poly2sym(Yin_dc,s);
        Yin = Yin_num/Yin_denum;
        %disp(Yin)
        %Yin = Ji2/Yin - Jri2/poly1;     %modify Yin here 
        % calculation of the Bi, check this step
        Bi = double(imag(subs(Yin, s, Tz(i+1))));      %error: division by zero
        %Bi = double(Bi);
        %disp(Bi);
        B = [B, Bi];    % append a new element
        Yin = Yin - 1i*Bi;
        count = count + 1;
        % disp(count)

    end
    %disp(Yin)
    %   extraction of the last dangling resonator, and the last Jrn
    %disp(Yin)
    [Yin_num,Yin_denum] = numden(Yin);
    Yin_num_c = double(sym2poly(Yin_num)); %numerator coefficients
    Yin_denum_c = double(sym2poly(Yin_denum)); %denumerator coefficients
    [res,poles,k] = residue(Yin_denum_c,Yin_num_c);
    % disp(res)
    % disp(poles)
    % disp(k)
    Jri2 = double(real(res(1)));    %the last section
    Jr = [Jr Jri2];
    Yin = k;    %check this part
    % disp(Yin)
    % debug the following part, check this
    % Bi = Ji2*imag(double(limit(Y_inverse,s,inf)));
    % J_last = Ji2*real(double(limit(Y_inverse,s,inf)))


    %assume the Gl = 1, calculate Bn and Bl
    %assumption making Bl = 0
    GL = 1/real(Yin);
    % disp(GL)
    %check the value of this
    BL_positive = sqrt((GL-(GL*GL)*real(Yin))/(real(Yin)));
    % disp(BL_positive)
    %disp(BL_positive)   % BL could be imaginary, could be solved by modify Bl
    BL_negative = -sqrt((GL-(GL*GL)*double(real(Yin)))/(double(real(Yin))));
    % disp(BL_negative)
    Bn1 = double(imag(Yin) + BL_positive/(BL_positive^2 + GL*GL));
    Bn2 = double(imag(Yin) + BL_negative/(BL_negative^2 + GL*GL));


end