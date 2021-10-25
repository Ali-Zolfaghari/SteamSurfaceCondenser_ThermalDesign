%***************************************************************************************************
%*   Calculate number of tube and length of them in steam surface condenser by presented code.
%*   I take no responsibilities for any errors in the code or damage thereby.
%*   Please notify me at zolfaghari1992iut@gmail.com if the code is used in any type of application.
%***************************************************************************************************
%*   Developer   : Ali Zolfaghari Sichani (20-06-2020)
%***************************************************************************************************
%*   References  : 
%* [1] “Standards for Steam Surface Condensers”, Heat Exchange Institute, 1995.
%* [2] “Problems Affecting Condenser Performance in Dr. Sherief Power Station”, said Ali, Basil Hassan, 2010, University of Khartoum.
%* [3] “The Computer Aided Design of Steam Surface Condensers, Jang”, JY, Leu, JS, Design and Operation of Heat Exchangers, 1992.
%* [4] “An Introduction to Chemical Cleaning of Industrial Water Systems”, Guyer, J Paul, 2018.
%* [5] “Heat exchangers: selection, rating, and thermal design”, Kakac, Sadik, Liu, Hongtan, Pramuanjaroenkij, Anchasa, 2020.
%* [6] “Compact heat exchangers”, Zohuri, Bahman, 2017.
%* [7] “Steam condensation analysis in a power plant condenser”, Drożyński, Zbigniew, Archives of thermodynamics, 2018.
%***************************************************************************************************
%*   Inputs      :
%*   DesignType  (1: HEI Standard - 2: Sadik Reference Book )
%*   other input parameters are clear by their names and in SI units.
%*   Outputs      :
%*   Npp_tube      (Number of Tube per Pass  )  (-)
%*   Lpp_tube      (Length of Tube per Pass  )  (m)
%***************************************************************************************************

    
clear,clc
close all
format long
format compact
load('LibData');

% =========================================================================
% INPUT %
% =========================================================================

DesignType = 2;


Mdot_steam = 5.0;
InletTemp_water = 22.0;
OutletTemp_water = 32.0;
InletTemp_steam = 60.0;
Velocity_water = 1.7;
Quality_steam = 100.0;
OverDesign_factor = 30.0;
Eta_pump = 85.0;

OD_tube = 0.0254;
BWG_tube = 16;
K_tube = 112.0;
Rfi_tube = 0.0002;
Rfo_tube = 0.0001;
Nc_tube = 25;
Type_tube = 30.0;
Np_tube = 2;
Pitch_tube = 0.03386;
Cleanliness_tube = 0.47;
Baffle_space = 0.08;

% =========================================================================
% SOLVER %
% =========================================================================

Gravity = 9.81;

fi = find( Fm(:,1) == K_tube );
fj = find( Fm(1,:) == BWG_tube );
fk = find( R2(1,:) == BWG_tube );
Thk_tube = Fm(2,fj)*0.0254;

BulkTemp_water = 0.5*(InletTemp_water+OutletTemp_water);
LatentEnthalpy_steam = 1000.0*(XSteam('hV_T',InletTemp_steam)-XSteam('hL_T',InletTemp_steam));
OutletEnthalpy_steam = XSteam('hL_T',InletTemp_steam)*1000.0;
InletEnthalpy_steam = OutletEnthalpy_steam+0.01*Quality_steam*LatentEnthalpy_steam;

Cp_water = XSteam('CpL_T',BulkTemp_water)*1000.0;
Rho_water = XSteam('rhoL_T',BulkTemp_water);
Mu_water = XSteam('my_PT',1,BulkTemp_water);
K_water = XSteam('tcL_T',BulkTemp_water);
Pr_water = (Mu_water*Cp_water)/K_water;

Rho_steam = XSteam('rhoL_T',InletTemp_steam);
Mu_steam = XSteam('my_PT',1,InletTemp_steam);
K_steam = XSteam('tcL_T',InletTemp_steam);

ID_tube = OD_tube-2.0*Thk_tube;
MD_tube = (OD_tube-ID_tube)/log(OD_tube/ID_tube);

LMTD = (OutletTemp_water-InletTemp_water)/log((InletTemp_steam-InletTemp_water)/(InletTemp_steam-OutletTemp_water));

Heat_condenser = Mdot_steam*(InletEnthalpy_steam-OutletEnthalpy_steam);

Mdot_water = Heat_condenser/(Cp_water*(OutletTemp_water-InletTemp_water));

Npp_tube = fix((4.0*Mdot_water)/(Velocity_water*Rho_water*pi*ID_tube*ID_tube));

if ( DesignType == 1)
    
    CFm = Fm(fi,fj);
    
    qx = zeros(size(U,1)-1,size(U,2)-1);
    qy = zeros(size(U,1)-1,size(U,2)-1);
    qv = zeros(size(U,1)-1,size(U,2)-1);
    for i = 2:size(U,1)
        for j = 2:size(U,2)
            qx(i-1,j-1) = U(1,j);
            qy(i-1,j-1) = U(i,1);
            qv(i-1,j-1) = U(i,j);
        end
    end
    
    if ( Np_tube == 1 )
        CREa = fitresult_RE1{1}(Velocity_water*3.28084);
        CREb = fitresult_RE1{2}(Velocity_water*3.28084);
        CREc = fitresult_RE1{3}(Velocity_water*3.28084);
    elseif ( Np_tube == 2 )
        CREa = fitresult_RE2{1}(Velocity_water*3.28084);
        CREb = fitresult_RE2{2}(Velocity_water*3.28084);
        CREc = fitresult_RE2{3}(Velocity_water*3.28084);
    elseif ( Np_tube == 3 )
        CREa = fitresult_RE3{1}(Velocity_water*3.28084);
        CREb = fitresult_RE3{2}(Velocity_water*3.28084);
        CREc = fitresult_RE3{3}(Velocity_water*3.28084);
    elseif ( Np_tube == 4 )
        CREa = fitresult_RE4{1}(Velocity_water*3.28084);
        CREb = fitresult_RE4{2}(Velocity_water*3.28084);
        CREc = fitresult_RE4{3}(Velocity_water*3.28084);
    end

    CRt = interp2(x_Rt,z_Rt,y_Rt,Velocity_water*3.28084,OD_tube/0.0254,'linear');
    
    CR2 = interp1(R2(2:end,1),R2(2:end,fk),OD_tube/0.0254,'linear');
    
    CR1 = interp1(R1(:,1),R1(:,2),1.8*BulkTemp_water+32.0,'linear');
    
    U1 = interp2(qx,qy,qv,Velocity_water*3.28084,OD_tube/0.0254,'linear');
    
    CFw = interp1(Fw(:,1),Fw(:,2),1.8*InletTemp_water+32.0,'linear');
    
    CFs = interp1(Fs(:,1),Fs(:,2),OD_tube/0.0254,'linear');
    
    A_tube = ((0.001*Heat_condenser*3412.14)/(U1*CFw*CFm*Cleanliness_tube*1.8*LMTD))*(1.0+0.01*OverDesign_factor);
    
    Lt_tube = (A_tube/(Npp_tube*CFs))/3.28084;
    
    Lpp_tube = Lt_tube/Np_tube;
    
    Dp_water = 2989.0669*(Lt_tube*3.28084*(CRt*CR2*CR1)+(CREa+CREb+CREc));
    
elseif ( DesignType == 2)
    
    Re_water = (Velocity_water*Rho_water*ID_tube)/Mu_water;
    
    fe_water = (1.58*log(Re_water)-3.28)^(-2.0);
    
    Nu_water = (0.5*fe_water*(Re_water-1000.0)*Pr_water)/(1.07+12.7*sqrt(0.5*fe_water)*((Pr_water^(2.0/3.0))-1.0));
    
    h_water = (Nu_water*K_water)/ID_tube;
    Rt = Rfo_tube+((1.0/h_water)+Rfi_tube)*(OD_tube/ID_tube)+(Thk_tube/K_tube)*(OD_tube/MD_tube);
    hs_steam = 0.728*(((Rho_steam*Rho_steam*Gravity*LatentEnthalpy_steam*K_steam*K_steam*K_steam)/(Mu_steam*OD_tube))^0.25)*(Nc_tube^(-1/6));
    
    Error = 1.0;
    DTw_1 = 10.05;
    while ( Error > 0.0001 )
        h_steam = hs_steam*(DTw_1^-0.25);
        U = 1.0/(Rt+(1.0/h_steam));
        DTw_2 = (InletTemp_steam-InletTemp_water)*(1.0-Rt*U);
        Error = abs(DTw_2-DTw_1);
        DTw_1 = DTw_2;
    end
    Ui = U;
    
    Error = 1.0;
    DTw_1 = 10.0;
    while ( Error > 0.0001 )
        h_steam = hs_steam*(DTw_1^-0.25);
        U = 1.0/(Rt+(1.0/h_steam));
        DTw_2 = (InletTemp_steam-OutletTemp_water)*(1.0-Rt*U);
        Error = abs(DTw_2-DTw_1);
        DTw_1 = DTw_2;
    end
    Uo = U;
    
    Um = 0.5*(Ui+Uo);
    
    A_tube = (Heat_condenser/(Um*LMTD))*(1.0+0.01*OverDesign_factor);
    
    Lt_tube = A_tube/(Npp_tube*pi*OD_tube);
    
    Lpp_tube = Lt_tube/Np_tube;
    
    CL_tube = 1.00;
    if ( Type_tube == 60.0 || Type_tube == 30.0 )
        CL_tube = 0.87;
    end
    
    CPT_tube = 0.93;
    if ( Np_tube == 2 )
        CPT_tube = 0.90;
    elseif ( Np_tube == 3 )
        CPT_tube = 0.85;
    end
    
    PR_tube = Pitch_tube/OD_tube;
    
    OD_shell = (0.637*sqrt(CL_tube/CPT_tube)*sqrt((A_tube*PR_tube*PR_tube*OD_tube)/(Lt_tube)))+2.0*Baffle_space;
    
    Dp_water = (4.0*fe_water*(Lt_tube/ID_tube)+4.0*Np_tube)*(0.5*Rho_water*Velocity_water*Velocity_water);
    
    P_pump = (Mdot_water*Dp_water)/(Rho_water*0.01*Eta_pump);
    
end 

    % =========================================================================
    % PRINT %
    % =========================================================================
    
    fprintf('=============================================================\n');
    fprintf('Number of Tube per Pass  :  %12d         \n',Npp_tube);
    fprintf('Length of Tube per Pass  :  %12.3f  (m)  \n',Lpp_tube);
    fprintf('Water Pressure Drop      :  %12.3f  (kPa)\n',Dp_water/1000.0);
    fprintf('=============================================================\n\n');
