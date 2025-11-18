fL=1930e6; fH=1995e6;   % lower band, upper band
N = 7;  % order of the filter
%Tz_tx = [1.44,-2.5,1.6,-1.65,1.6,-2.5,1.74];
Tz_tx = [2.26,-2.26,1.63,-1.42,1.63,-2.26,2.26];
%Tz_rx = [2.632,-2.223,2.079,-2.074,2.080,-2.228,2.599];
Tz_rx = [-1.7,1.97,-2.5,3,-3.3,4,-1.2];
% Tz1 = [2.26,-2.26,1.63,-1.61,1.63,-2.26,2.26];
RL = 15;    % return loss
xx=-3:0.01:10;  % low pass plting range
f0 = sqrt(fL*fH);
FBW = (fH-fL)/f0;
alpha = f0/(fH-fL);


%   step 1: define the filter function
[S11_Zeros,E_Zeros,LP_w,e,er] = Define_LP_filter(fL,fH,Tz_tx,N,RL,xx);
%   LP_w is the predefined Tz
%   S11_z is the zeros of F(s)
%   E_z is the zeros of E(s)

%   Construct F(s),P(s),E(s)
syms s
F = 1;  % initial F
P = 1;  % initial P
E = 1;  % initial E

%   construct P
for i = 1:length(LP_w)
    P = P*(s-LP_w(i));
end
%   construct F
for j = 1:length(S11_Zeros)
    F = F*(s-S11_Zeros(j));
end
%   construct E
for k = 1:length(E_Zeros)
    E = E*(s- E_Zeros(k));
end


%   Coefficients of F,P,E
Pc = sym2poly(P);
Fc = sym2poly(F);
Ec = sym2poly(E);
%   Construct S11 and S12, input impedance
S11 = F/(E*er);
S12 = P/(E*e);
[S11_num,S11_denum] = numden(S11);
[S12_num,S12_denum] = numden(S12);
zeros_sym = solve(S12_num == 0, s);


%   step 2: calculate the input phase
input_phase = double((subs(S11_denum,s,LP_w(1)))/(subs(S11_num,s,LP_w(1))));    
S11 = S11*input_phase;
Y_in = (1-S11)/(1+S11); 
[Y_in_num,Y_in_denum] = numden(Y_in);
Y_in =  Y_in_num/Y_in_denum;
Yin_num_c = double(sym2poly(Y_in_num)); %numerator coefficients
Yin_denum_c = double(sym2poly(Y_in_denum)); %denumerator coefficients
[res,poles] = residue(Yin_denum_c,Yin_num_c);   %seems ok for this step

%   step 3: start the parameter extraction
Ji2 = 1;    %design parameter
[B,Jr,Bn1,Bn2] = Parameter_extraction(Y_in,Ji2,LP_w,1);
B = real([B, Bn1]);
%   debug check the order of B
%   se-sh-se-sh-....
%   the order of B is same with the order of Jr

%   step 4: do the LP-BP transformation
b = -imag(LP_w);
w0 = 2*pi*f0;
%   debug this function，the output is negative sometimes
[Lse_a,Cse_a,Cse_0,Lsh_a,Csh_a,Csh_0] = LP_BP_transformation(B,Jr,b,50,alpha,w0,N);







