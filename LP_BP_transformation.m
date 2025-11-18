%   Transformation from LP to BP

function  [Lse_a,Cse_a,Cse_0,Lsh_a,Csh_a,Csh_0] = LP_BP_transformation(B,Jr,Tz,Z0,alpha,wo,N)
%   Output
%   LA,CA,CO is the lumped element value start from the first resonator to
%   the last one
%   Z0:  characteristic impedance
%   alpha: FBW
%   wo: center frequency
%   N: order of the filter




%   start from series
%   Split the NRN into series or shunt
    B_se = [];
    B_sh = [];

    B_se = B(1:2:end);   % 奇数位置
    B_sh = B(2:2:end);   % 偶数位置
    % for i = 1:length(N)
    %     if (rem(i,2) == 1)
    %         B_se = [B_se, B(i)];
    %     else
    %         B_sh = [B_sh,B(i)];
    %     end
    % 
    % end
     % disp(B_sh)
     % disp(B_se)

%   Split the Jr into series or shunt
    Jr_se = [];
    Jr_sh = [];
    Jr_se = Jr(1:2:end);   % 奇数位置
    Jr_sh = Jr(2:2:end);   % 偶数位置
    % for i = 1:length(N)
    %     if (rem(i,2) == 1)
    %         Jr_se = [Jr_se, Jr(i)];
    %     else
    %         Jr_sh = [Jr_sh,Jr(i)];
    %     end
    % 
    % end 
     % disp(Jr_sh)
     % disp(Jr_se)

%   Split the Tz into series or shunt
    Tz_se = [];
    Tz_sh = [];
    Tz_se = Tz(1:2:end);   % 奇数位置
    Tz_sh = Tz(2:2:end);   % 偶数位置
    % for i = 1:length(N)
    %     if (rem(i,2) == 1)
    %         Tz_se = [Tz_se, Tz(i)];
    %     else
    %         Tz_sh = [Tz_sh,Tz(i)];
    %     end
    % end

     % disp(Tz_sh)
     % disp(Tz_se)




    %shunt
    Xsh_0 = []; %list of LP Xsh_0
    Xsh_m = []; %list of LP Xsh_m 
    Lsh_m = []; %list of LP Lsh_m
    for i = 1:length(B_sh)
        Xsh0 = -1/B_sh(i);
        Xsh_0 = [Xsh_0,Xsh0];
    end
    
    for j = 1:length(B_sh)
        Xshm = Tz_sh(j)/Jr_sh(j);
        Xsh_m = [Xsh_m,Xshm];
    end
    
    for k = 1:length(B_sh)
        Lshm  = 1/Jr_sh(k);
        Lsh_m = [Lsh_m,Lshm];
    end

     % disp(Lsh_m);
     % disp(Xsh_m);
     % disp(Xsh_0);

    
    %series
    Xse_0 = []; %list of LP Xsh_0
    Xse_m = []; %list of LP Xsh_m 
    Lse_m = []; %list of LP Lsh_m
    % assume Jm2 = 1
    for i = 1:length(B_se)
        Xse0 = B_se(i);
        Xse_0 = [Xse_0,Xse0];
    end
    
    for j = 1:length(B_se)
        Lsem = (B_se(j))^2/Jr_se(j);
        Lse_m = [Lse_m,Lsem];
    end
    
    for k = 1:length(B_se)
        %disp(Tz_se);
        inter = Tz_se(k)*B_se(k)/Jr_se(k);
        Xsem  = B_se(k)*(inter-1);
        Xse_m = [Xse_m,Xsem];
    end

    % disp(Lse_m);    %debug this
    % disp(Xse_m);
    % disp(Xse_0);
    

    %   transform to BP equivalent circuit
    %   series 
    Lse_a = [];
    Cse_a = [];
    Cse_0 = [];

    for i = 1:length(Lse_m)
        Lsea = 0.5*(2*alpha*Lse_m(i) + Xse_m(i))*Z0/wo;
        Csea = 2/(wo*Z0*(2*alpha*Lse_m(i) - Xse_m(i)));
        Cse0 = -1/(wo*Z0*Xse_0(i));

        Lse_a = [Lse_a,Lsea];
        Cse_a = [Cse_a,Csea];
        Cse_0 = [Cse_0,Cse0];
    end

    %   shunt
    Lsh_a = [];
    Csh_a = [];
    Csh_0 = [];

    for i = 1:length(Lsh_m)
        Lsha = 0.5*(2*alpha*Lsh_m(i) + Xsh_m(i))*Z0/wo;
        Csha = 2/(wo*Z0*(2*alpha*Lsh_m(i) - Xsh_m(i)));
        Csh0 = -1/(wo*Z0*Xsh_0(i));

        Lsh_a = [Lsh_a,Lsha];
        Csh_a = [Csh_a,Csha];
        Csh_0 = [Csh_0,Csh0];
    end
    

end