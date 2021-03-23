%AnisotropicConduction
%Calculation of conduction through structure of coupled traxels

close all
clear all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% In this script we use (Ntraxels x 2) matrices to define inputs and
% outputs which establish the boundary conditions. The meaning of the
% various values is:
% First element: integer [1,2,3], select kind of boundary
% 1=voltage, 2=current=0, 3=connection to other traxel
% Second element: real 
% case 1=voltage value
% case 2=redundant
% case 3=other traxel number
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

RunModel;


% Parameters
[Ntraxels,Width,Height,Length,Rho,Sigma,Uin,Uou]=InitParameters;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Initialize calculations to represent Transducers Paper structure 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Ntraxels,Width,Height,Length,Rho,Sigma,Uin,Uou]=InitParameters
    %Speficy geometric parameters
    Ntraxels = 19;                 %Number of traxels in the structure  
    Width = 0.8E-3 ;               %Width of traxels in m
    Height = 200e-6;               %Height of traxels in m
    Length = 15E-3;                %Length of traxels in m
    Area = Width*Height;           %Crossectional area of traxel in m^2

    %Specify material parameters
    Rho = 2.8;                    %Volume Resistivity of filament in Ohm*m
    Sigma = 2e-3;                 %Surface resistivity in Ohm*m^2
    
    % Anisotropy Ratio
    Ani = Width*Rho/(Width*Rho + Sigma);
    assignin('base','Ani',Ani);
    
    % Aspect Ratio
    AR = Ntraxels*Width/Length;
    assignin('base','AR',AR);
    
    Uhig=1;                     %High voltage
    Ulow=0;                      %Low voltage


    %Define connections
    %Specify inputs/outputs
    %First element: integer [1,2,3], select kind of boundary
    %1=voltage, 2=current=0, 3=connection to other traxel
    %Second element: real 
    %case 1=voltage value
    %case 2=redundant
    %case 3=other traxel number
    

    % Opposite Side Leads
    [Uin,Uou]=AllOpenConnections(Ntraxels);
    [Uin]=Meander(Uin,2,Ntraxels);
    [Uou]=Meander(Uou,1,Ntraxels);
    Uin( 1: 1,1:2)=[1 Uhig];
    Uou( end:end,1:2)=[1 Ulow];
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to calculate and simulated data measured by Alexander Dijkshoorn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function RunModel
    %Specify geometric and material parameters and boundary conditions
    
    [N,W,H,L,Rho,Sigma,Uin,Uou]=InitParameters;
    
    %Some plot parameters
    Npoints=200;
    %Initialize Voltage matrix
    U(1:N,1:Npoints)=0;
    %Initialize Current matrix
    I(1:N,1:Npoints)=0;
    %Initialize Power matrix
    P(1:N,1:Npoints)=0;
    %Initialize z-positions
    z(1:Npoints)=10;
    
    
    [z,U,I,Jx,P,Pdens,Z]=CalcAllProfiles(N,L,W,H,Rho,Sigma,Uin,Uou,Npoints);
    PlotAll(W,z,real(U),real(I),Jx,P,Pdens,Z);
    assignin('base', 'position', z);
    assignin('base', 'R', Z);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to have all input and output disconnected
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Uin,Uou]=AllOpenConnections(N)
    Uin(1:N,1:1)=2;
    Uin(1:N,2:2)=0;
    Uou(1:N,1:1)=2;
    Uou(1:N,2:2)=0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to have all input and output connected like meander
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Uin,Uou]=AllMeanderInput(N)
    P=floor(N/2);
    for m=1:P
        Uin(2*m:2*m    ,1:2)=[3 2*m+1];
        Uin(2*m+1:2*m+1,1:2)=[3 2*m];
    end
    for m=1:P
        Uou(2*m-1:2*m-1,1:2)=[3 2*m];
        Uou(2*m:2*m    ,1:2)=[3 2*m-1];
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to have defined inputs or outputs connected
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Uout]=Meander(Uinp,First,Last)
    Uout=Uinp;
    for m=First:Last
        dm=1-2*round(mod(m-First,2)/2);
        Uout(m:m,1:2)=[3 m+dm];
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to have all input connected to Uhig
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Uin]=AllInputConnected(N,Uhig)
    Uin(1:N,1:1)=1;
    Uin(1:N,2:2)=Uhig;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to have all output connected to Ulow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Uou]=AllOutputConnected(N,Ulow)
    Uou(1:N,1:1)=1;
    Uou(1:N,2:2)=Ulow;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function initialize some variables from the physical input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Area,Zint,Rlat,Gamma]=Initialize(Width,Height,Rho,Sigma)
    %First calculate some needed help variables
    Area=Width*Height;                   %Crossectional area of traxel in m^2
    Rint=Sigma/Height;                   %Interface resistance in Ohm*m
    Zint=Rint;                           %Interface impedance in Ohm/m
    Rlat=Rho/Area;                       %Lateral resistance in Ohm/m
    
    Gamma = Rho*Width./(Rho*Width+Sigma);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to calculate the eigenvalues (eVac), eigenvectors (eVec) 
% coefficientss (Alpha) and impedance (Z) for s structure of Ntrax
% structure length (L), number of traxels (N), gamma value (G)
% and boundary conditions (Uin, Uou)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Lambda,eVec,Alpha,Z]=CalcSolution(Ntrax,Length,Width,Height,Rho,Sigma,Uin,Uou)

    Nlambda=2*Ntrax;
    B(1:Nlambda,1:Nlambda)=0;
    X(1:Nlambda)=0;
    L2=0.5*Length;

    [Area,Zint,Rlat,Gamma]=Initialize(Width,Height,Rho,Sigma);

    %Fill matrix A with the coefficients from the coupled equations
    A(1:Ntrax,1:Ntrax)=0;
    for m=1:Ntrax
        for k=1:Ntrax
            if m==k 
                if (m>1) && (m<Ntrax) 
                    A(m,k)=-2*Gamma/(Width^2);
                else
                    A(m,k)=-1*Gamma/(Width^2);
                end
            else
                if abs(m-k)==1
                    A(m,k)=Gamma/(Width^2);
                end
            end
        end
    end

    %Calculate eigenvalues and eigenvectors. Note: it seems eig(A) calculates
    %the eigenvalues of (A-lambda*I).
    [eVec,eVal]=eig(-A);

    %Now sort the indices and eigenvectors
    [d,ind]=sort(diag(eVal));
    eVal=eVal(ind,ind);
    eVec=eVec(:,ind);

    % Generate the Lambda's from the eigenvalues.
    for m=1:Ntrax
        Lambda(2*m-1)= sqrt(eVal(m,m));
        Lambda(2*m  )=-sqrt(eVal(m,m));
    end

    % Note: m(i,j) means m(row=i, column=j)
    % Input conditions for x=0
    for k=1:Ntrax           %Go over all layers
        switch Uin(k:k,1:1)
            case 1 %Voltage specified
                X(k)=Uin(k:k,2:2);
                B(k,1)= eVec(k,1);
                B(k,2)=-eVec(k,1)*L2;
                for m=3:Nlambda B(k,m)=eVec(k,ceil(m/2))*exp(-Lambda(m)*L2);
                end
            case 2 %Current is 0 condition
                X(k)=0;
                B(k,1)=0;
                B(k,2)=eVec(k,1);
                for m=3:Nlambda B(k,m)=eVec(k,ceil(m/2))*Lambda(m)*exp(-Lambda(m)*L2);
                end
            case 3 %Traxel inputs are connected
                X(k)=0;
                cnc=Uin(k:k,2:2);   %cnc=connecting channel
                if cnc>k            %apply voltage condition U_cnc-U_k=0. 
                    B(k,1)= (eVec(cnc,1)-eVec(k,1));
                    B(k,2)=-(eVec(cnc,1)-eVec(k,1))*L2;
                    for m=3:Nlambda B(k,m)=(eVec(cnc,ceil(m/2))-eVec(k,ceil(m/2)))*exp(-Lambda(m)*L2);
                    end 
                end
                if cnc<k      %Apply current condition I_k+I_cnc=0
                    B(k,1)=0;
                    B(k,2)=eVec(cnc,1)+eVec(k,1);
                    for m=3:Nlambda B(k,m)=(eVec(cnc,ceil(m/2))+eVec(k,ceil(m/2)))*Lambda(m)*exp(-Lambda(m)*L2);
                    end
                end
        end 
    end


    % Note: m(i,j) means m(row=i, column=j)
    % Output conditions for x=L
    for k=1:Ntrax
        s=k+Ntrax;
        switch Uou(k:k,1:1)
            case 1 %Input voltage specified, current <> 0
                X(s)=Uou(k:k,2:2);
                B(s,1)=eVec(k,1);
                B(s,2)=eVec(k,1)*L2;
                for m=3:Nlambda 
                    B(s,m)=eVec(k,ceil(m/2))*exp(Lambda(m)*L2);
                end
            case 2 %No current, output not connected
                X(s)=0;
                B(s,1)=0;
                B(s,2)=eVec(k,1);
                for m=3:Nlambda 
                    B(s,m)=eVec(k,ceil(m/2))*Lambda(m)*exp(Lambda(m)*L2);
                end
            case 3 %Traxel inputs are connected
                X(s)=0;
                cnc=Uou(k:k,2:2);   %cnc=connecting channel
                if cnc>k            %apply voltage condition U_k-U_cnc=0. 
                    B(s,1)=(eVec(cnc,1)-eVec(k,1));
                    B(s,2)=(eVec(cnc,1)-eVec(k,1))*L2;
                    for m=3:Nlambda 
                        B(s,m)=(eVec(cnc,ceil(m/2))-eVec(k,ceil(m/2)))*exp(Lambda(m)*L2);
                    end 
                end
                if cnc<k      %Apply current condition I_k+I_cnc=0
                    B(s,1)=0;
                    B(s,2)=eVec(cnc,1)+eVec(k,1);
                    for m=3:Nlambda 
                        B(s,m)=(eVec(cnc,ceil(m/2))+eVec(k,ceil(m/2)))*Lambda(m)*exp(Lambda(m)*L2);
                    end
                end

        end
    end

    Y=transpose(X);
    Alpha=linsolve(B,Y);
    assignin('base','Y',Y)
    assignin('base','B',B)

    assignin('base', 'Alpha', Alpha);
    assignin('base', 'eVal', eVal);
    assignin('base', 'eVec', eVec);

    %Calculate impedance. First the current at the input and output
    for k=1:Ntrax           %k runs over the number of traxels
        Iin(k)=Alpha(2)*eVec(k,1);
        Iou(k)=Alpha(2)*eVec(k,1);
        for q=3:Nlambda 
            Iin(k)=Iin(k)+Alpha(q)*eVec(k,ceil(q/2))*Lambda(q)*exp(-Lambda(q)*L2);
            Iou(k)=Iou(k)+Alpha(q)*eVec(k,ceil(q/2))*Lambda(q)*exp( Lambda(q)*L2);
        end
    end

    Iin=-(Area/Rho)*Iin;
    Iou=-(Area/Rho)*Iou;

    Uhig=max(X);
    Ulow=min(X);
    Iinp=max(max(Iin),sum(Iin));
    if real(Iinp)<0
        Iinp=-1*Iinp;
    end
    Z=(Uhig-Ulow)/Iinp;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to calculate all voltages (U), currents (I) and powers (P) from
% the input structure length (L), number of traxels (N), gamma value (G)
% and boundary conditions (Uin, Uou)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [z,U,I,Jx,P,Pdens,Z]=CalcAllProfiles(Ntrax,Length,Width,Height,Rho,Sigma,Uin,Uou,Npoints)

    L2=0.5*Length;
    for m=1:Npoints 
        z(m)=(m-1)*Length/(Npoints-1); 
    end
    P(1:Ntrax,1:Npoints)=0;
    deltaZ=Length/Npoints;
    deltaA=deltaZ*Height;
    Nlambda=2*Ntrax;

    [Area,Zint,Rlat,Gamma]=Initialize(Width,Height,Rho,Sigma);

    [Lambda,eVec,Alpha,Z]=CalcSolution(Ntrax,Length,Width,Height,Rho,Sigma,Uin,Uou);

    % Calculate the voltages in the traxels as function of position
    for k=1:Ntrax           %k runs over the number of traxels
        for m=1:Npoints     %m runs over the position in each traxel
            zp=z(m)-L2;
            U(k,m)=Alpha(1)*eVec(k,1)+Alpha(2)*eVec(k,1)*zp;
            for q=3:Nlambda
                U(k,m)=U(k,m)+Alpha(q)*eVec(k,ceil(q/2))*exp(Lambda(q)*zp);
            end
        end
    end

    assignin('base','U',U)       % write U voltage to workspace

    % Calculate the currents in the traxels as function of position
    for k=1:Ntrax           %k runs over the number of traxels
        for m=1:Npoints     %m runs over the position in each traxel
            zp=z(m)-L2;
            I(k,m)=Alpha(2)*eVec(k,1);
            for q=3:Nlambda
                I(k,m)=I(k,m)+Alpha(q)*eVec(k,ceil(q/2))*Lambda(q)*exp(Lambda(q)*zp);
            end
        end
    end

    I=-(Area/Rho)*I;
    Jx = I/(Area);
    assignin('base', 'Jx', I/Area);


    % Calculate power-density in traxels as function of position
    dZc = 1/(Height*Length/Npoints)*( Rho*Width+Sigma) / (1);   % inter-traxel impedance

    for k=1:Ntrax           %k runs over the number of traxels
        for m=1:Npoints     %m runs over the position in each traxel
            if m>1
                P(k,m)= abs((U(k,m)-U(k,m-1))*I(k,m));            
            elseif m<=1
                P(k,m)= abs((U(k,m+1)-U(k,m))*I(k,m));             
            end

            if k>1          
                P(k,m)=P(k,m)+1/2*abs((U(k,m)-U(k-1,m)).*(U(k,m)-U(k-1,m))./(dZc));          
            end
            if k<Ntrax
                P(k,m)=P(k,m)+1/2*abs((U(k,m)-U(k+1,m)).*(U(k,m)-U(k+1,m))./(dZc));  
            end
        end
    end

    assignin('base','P',P)
    Pdens = P/(Width*Height*Length/Npoints);        % [W/m^3] Power density at every position
    Ptot = sum(sum(P));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to do all the plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PlotAll(W,z,U,I,Jx,P,Pdens,Z)

    figure
    s1=['U_1(z)'];
    s2=['U_2(z)'];
    s3=['U_3(z)'];
    s4=['I_1(z)'];
    s5=['I_2(z)'];
    s6=['I_3(z)'];
    s7=['I_{tot}(z)'];

    col=['r' 'g' 'b' 'y' 'm' 'c' 'r--' 'g--' 'b--' 'y--' 'm--' 'c--'];
    [cr,cl]=size(col);

    lw=2; 
    FS=17;
    AFS=13;

    [Ntraxels,Npoints]=size(U);
    L=max(z);
    H=Ntraxels*W;
    shelp='F=%.0f Hz';
    sf=compose(shelp);
    shelp='Z=%.0f+j%.0f';
    sz=strcat(compose(shelp,real(Z),imag(Z)),' \Omega');

    subplot(2,3,1)
    hold on;
    grid on;
    for m=1:Ntraxels 
       plot(z,U(m:m,1:Npoints),col(mod(m-1,cl)+1),'linewidth',lw); 
    end

    title(strcat('Voltage vs position. '))
    xlabel('x [m]','fontsize',FS)
    ylabel('U_j(x) [V]','fontsize',FS);
    set(gca,'FontSize',AFS);
    ylim([0 1])
    xlim([min(z) max(z)])

    [Ntraxels,Npoints]=size(I);
    assignin('base','I',I);

    subplot(2,3,2)
    hold on;
    grid on;
    for m=1:Ntraxels
       plot(z,Jx(m:m,1:Npoints),col(mod(m-1,cl)+1),'linewidth',lw); 
    end

    title(strcat('Current Density x-direction vs position. '));
    xlabel('x [m]','fontsize',FS)
    ylabel('J_{x,j}(x) [A/m^2]','fontsize',FS);
    set(gca,'FontSize',AFS);
    xlim([min(z) max(z)])

    [Ntraxels,Npoints]=size(P);

    subplot(2,3,3)
    for m=1:Ntraxels
        semilogy(z,(Pdens(m:m,1:Npoints)),col(mod(m-1,cl)+1),'linewidth',lw);
       hold on;
        grid on;
    end

    title('Power Density vs position')
    xlabel('x [m]','fontsize',FS)
    ylabel('Q_j(x) [W/m^3]','fontsize',FS);
    set(gca,'FontSize',AFS);
    xlim([min(z) max(z)])

    MaxCol=2048;
    cm=jet(MaxCol); %Set color scales


    % Establish real-world dimensions
    RI=imref2d(size(U));
    RI.XWorldLimits=[0 L];
    RI.YWorldLimits=[0 H];
    Umin=min(U(:));
    Umax=max(U(:));
    x=W/2:W:H;

    subplot(2,3,4)
    hold on;
    grid on;
    imshow(real(U),RI,[Umin Umax]);
    colormap(gca, cm);

    colorbar;
    title('Voltage [V]');
    xlabel('x [m]','fontsize',FS);
    ylabel('y [m]','fontsize',FS);
    set(gca,'FontSize',AFS);

    % Establish real-world dimensions
    RI=imref2d(size(I));
    RI.XWorldLimits=[0 L];
    RI.YWorldLimits=[0 H];
    Jmax=max(Jx(:));
    Jmin=min(Jx(:));

    stitle=strcat('Current Density x-direction [A/m^2]');
    subplot(2,3,5)
    hold on;
    grid on;
    imshow(real(Jx),RI,[Jmin Jmax]);
    colormap(gca, cm);
    colorbar;
    title(stitle);
    xlabel('x [m]','fontsize',FS);
    ylabel('y [m]','fontsize',FS);
    set(gca,'FontSize',AFS);
    set(gca,'FontSize',AFS);


    % Establish real-world dimensions
    RI=imref2d(size(P));
    RI.XWorldLimits=[0 L];
    RI.YWorldLimits=[0 H];


    Pmaxd=max(Pdens(:));
    Pmind=min(Pdens(:));
    cm=hot(MaxCol);
    
    
    subplot(2,3,6)
    hold on;
    grid on;
    imshow((Pdens),RI,[(Pmind) (Pmaxd)]);
    colormap(gca, cm);
    colorbar;
    title('Power Density [W/m^3]');
    xlabel('x [m]','fontsize',FS);
    ylabel('y [m]','fontsize',FS);
    set(gca,'FontSize',AFS);
    set(gca,'FontSize',AFS);
    set(gca,'ColorScale','log')

end