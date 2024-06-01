% Matteo Parente, 20/07/2016, DIFFER Eindhoven
% Andrea Baldi, 30/05/2024 revised, Vrije Universiteit Amsterdam

% This script reads the relative permittivity data of two materials and
% uses Mie theory to compute the scattering, extinction, and absorption
% cross-sections for a spherical core-shell particle embedded in a lossless 
% medium. The permittivity files need be tab-delimited text files with
% three columns: energy (in eV), epsilon1, epsilon2.

clear all
close all

% Define default values
default_r_in = '25'; % Default core radius in nm
default_t_shell = '25'; % Default shell thickness in nm
default_ref_med = '1.333'; % Default refractive index of the surrounding medium (air = 1, water = 1.333)
default_E_min = '1'; % Default minimum energy in eV
default_E_max = '5'; % Default maximum energy in eV
default_Eval = '401'; % Default number of energy values to compute
default_nbes = '10'; % Default maximum order of the Bessel function

% Manual input
prompt = {'Core radius (nm)','Shell thickness (nm)','Refractive index of the surrounding medium','Minimum energy value','Maximum energy value','Number of energy values to compute','Maximum order of the Bessel function'};
dlg_title = 'Input parameters';
num_lines = 1;
def = {default_r_in, default_t_shell, default_ref_med, default_E_min, default_E_max, default_Eval, default_nbes};
answer = inputdlg(prompt,dlg_title,num_lines,def);
r_in=str2double(answer{1}); % Core radius in nm
t_shell=str2double(answer{2}); % Shell thickness in nm
r_out=r_in+t_shell; % Radius of the outer sphere in nm 
ref_med=str2double(answer{3}); % Refractive index of the surrounding medium
Emin=str2double(answer{4}); % Minimum energy in eV
Emax=str2double(answer{5}); % Maximum energy in eV
Eval=str2double(answer{6}); % Number of energy values to compute
nbes=str2double(answer{7}); % Maximum order of the Bessel function

% Fundamental physical constants in SI units
e = 1.60217646e-19;     % Elementary charge 
h = 6.626068e-34;       % h (Planck's constant)
hbar = 1.05457148e-34;  % hbar (Planck's constant divided by 2*Greek Pi)
me = 9.10938215e-31;    % Rest mass of the electron 
c = 2.99792458e8;       % Speed of light

% Select permittivities
% Save results
hMsgBox = msgbox('Select the permittivity of the core', 'modal');
uiwait(hMsgBox);
read_in=dlmread(uigetfile('*','Select a dielectric function file for the inner sphere')); %Loads a file in which are specified enertgy (eV), eps1 (real) and eps2 (imaginary) for the inner sphere
eV_in=read_in(:,1);
e1_in=read_in(:,2);
e2_in=read_in(:,3);
hMsgBox = msgbox('Select the permittivity of the shell', 'modal');
uiwait(hMsgBox);
read_out=dlmread(uigetfile('*','Select a dielectric function file for the outer sphere')); %Loads a file in which are specified enertgy (eV), eps1 (real) and eps2 (imaginary) for the outer sphere
eV_out=read_out(:,1);
e1_out=read_out(:,2);
e2_out=read_out(:,3);

% Calculations

%Interpolation of the data from the uploaded file in order to have a smooth data set
Energy=(Emin:(Emax-Emin)/(Eval-1):Emax)';
e1_in=interp1(eV_in,e1_in,Energy, 'spline');
e2_in=interp1(eV_in,e2_in,Energy, 'spline');
e1_out=interp1(eV_out,e1_out,Energy, 'spline');
e2_out=interp1(eV_out,e2_out,Energy, 'spline');

%Calculation of the total permittivity values for the inner and outer sphere
n_in=(((e1_in.^2 + e2_in.^2).^(1/2) + e1_in)./2).^(1/2);
k_in=(((e1_in.^2 + e2_in.^2).^(1/2) - e1_in)./2).^(1/2);
n_tilde_in=n_in+1i*k_in;
n_out=(((e1_out.^2 + e2_out.^2).^(1/2) + e1_out)./2).^(1/2);
k_out=(((e1_out.^2 + e2_out.^2).^(1/2) - e1_out)./2).^(1/2);
n_tilde_out=n_out+1i*k_out;
m_in=n_tilde_in./ref_med;
m_out=n_tilde_out./ref_med;

%Wavelength, wavenumber k and conversion of radius in meters
lambda=h*c./(e*Energy);
k=2*pi*ref_med./lambda;
radius_in=r_in*1e-9;
radius_out=r_out*1e-9;

%Size parameters and inizialization of the scattering elements
x=k.*radius_in;
y=k.*radius_out;
z1=m_in.*x;
z2=m_out.*x;
z2y=m_out.*y;
costx=sqrt(pi./(2*x));
costy=sqrt(pi./(2*y));
costz1=sqrt(pi./(2*z1));
costz2=sqrt(pi./(2*z2));
costz2y=sqrt(pi./(2*z2y));
sca_elem=0;
ext_elem=0;

%Beginning of the for cycle to calculate the scattering and extinction elements with iteration of the Bessel function
for n=1:nbes

% Bessel functions
jnz1=costz1.*(besselj(n+0.5,z1));
jnz2=costz2.*(besselj(n+0.5,z2));
jnz2y=costz2y.*(besselj(n+0.5,z2y));
jny=costy.*(besselj(n+0.5,y));
ynz1=costz1.*bessely(n+0.5,z1);
ynz2=costz2.*bessely(n+0.5,z2);
ynz2y=costz2y.*bessely(n+0.5,z2y);
yny=costy.*bessely(n+0.5,y);
hny=jny+1i.*yny;
psiz1=z1.*jnz1;
psiz2=z2.*jnz2;
psiz2y=z2y.*jnz2y;
psiy=y.*jny;
chiz1=-z1.*ynz1;
chiz2=-z2.*ynz2;
chiz2y=-z2y.*ynz2y;
xiy=y.*hny;

%Derivatives
psiz1der=z1.*costz1.*besselj(n-0.5,z1)-n.*jnz1;
psiz2der=z2.*costz2.*besselj(n-0.5,z2)-n.*jnz2;
psiz2yder=z2y.*costz2y.*besselj(n-0.5,z2y)-n.*jnz2y;
psiyder=y.*costy.*besselj(n-0.5,y)-n.*jny;
chiz1der=-z1.*costz1.*bessely(n-0.5,z1)+n.*ynz1;
chiz2der=-z2.*costz2.*bessely(n-0.5,z2)+n.*ynz2;
chiz2yder=-z2y.*costz2y.*bessely(n-0.5,z2y)+n.*ynz2y;
xiyder=0.5.*hny+(0.5.*y.*costy.*(besselj(n-0.5,y)+1i.*bessely(n-0.5,y)-besselj(n+1.5,y)-1i.*bessely(n+1.5,y)));

%Coefficients
An=(m_out.*psiz2.*psiz1der-m_in.*psiz2der.*psiz1)./(m_out.*chiz2.*psiz1der-m_in.*chiz2der.*psiz1);
Bn=(m_out.*psiz1.*psiz2der-m_in.*psiz2.*psiz1der)./(m_out.*chiz2der.*psiz1-m_in.*psiz1der.*chiz2);
an=(psiy.*(psiz2yder-An.*chiz2yder)-m_out.*psiyder.*(psiz2y-An.*chiz2y))./(xiy.*(psiz2yder-An.*chiz2yder)-m_out.*xiyder.*(psiz2y-An.*chiz2y));
bn=(m_out.*psiy.*(psiz2yder-Bn.*chiz2yder)-psiyder.*(psiz2y-Bn.*chiz2y))./(m_out.*xiy.*(psiz2yder-Bn.*chiz2yder)-xiyder.*(psiz2y-Bn.*chiz2y));
sca_elem=sca_elem+((2*n+1).*(abs(an).^2+abs(bn).^2));
ext_elem=ext_elem+((2*n+1).*real(an+bn));
end

% Scattering, extinciton, and absorption cross sections and efficiencies
Csca = 2*pi./(k.^2).*sca_elem;
Cext = 2*pi./(k.^2).*ext_elem;
Cabs = Cext-Csca;
Qsca = Csca./(pi*radius_out^2);
Qext = Cext./(pi*radius_out^2);
Qabs = Cabs./(pi*radius_out^2);

% Plot results
fontSize = 12;
figure(1);
plot(Energy, Cext, 'k', 'LineWidth', 1);
hold on;
plot(Energy, Csca, '--k', 'LineWidth', 1);
plot(Energy, Cabs, '-.k', 'LineWidth', 1);
xlabel('Energy (eV)', 'FontSize', fontSize);
ylabel('Cross sections (m^2)', 'FontSize', fontSize);
lgd = legend('Extinction', 'Scattering', 'Absorption', 'Location', 'NorthEast');
set(lgd, 'FontSize', fontSize, 'Box', 'off');
set(gca, 'FontSize', fontSize);
