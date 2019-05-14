%Writing up the pipeline for external FTMS users
%KLongnecker 2/25/09; edited EKujawinski 2/27/09; KLongnecker 9/2/09
%KLongnecker 8/9/10; KLongnecker 1/24/2013; KLongnecker 7/5/2013; KL
%4/28/2105 to correct error in findformula
%KL 11/29/2016 updating to use new database with corrected P valence 
%and SENIOR rule
%
% NOTE: The database required to run this code exceeds GitHub's limits on
% file sizes. The database is available by anonymous FTP at this location:
% ftp://ftp.whoi.edu/pub/users/klongnecker/LabFiles
% The file is also available from here:
% https://drive.google.com/a/whoi.edu/file/d/1pgpMlcg63C7pZO5bWPmyDbOGPvvYC59L/view?usp=sharing
% The name of the file is: 
%LongneckerKujawinski_fullCompoundList.2016.11.21.mat
%
%%contact information:
% Krista Longnecker
% Department of Marine Chemistry and Geochemistry
% Woods Hole Oceanographic Institution
% klongnecker@whoi.edu

%In the following documentation, anything which is not commented out 
%(lines beginning with '%') is something to be entered in a Matlab command
%window

%1. Align the peaks within the desired margin of error. For the WHOI
%instrument, we have been using 1 ppm. In the line below, the first '1'
%indicates that the code should use a 1 ppm margin of error for the peaks.
%The second '1' instructs Matlab to show you its progress - that can be
%changed to '0' if you don't want to have any information printed out in
%the command window.
%Note that code expects 'data' to be a matrix with multiple cells. There
%will be one cell per sample. Within each cell, there is a two column
%matrix. The first column is the list of m/z values in that sample and the
%second column is the peak heights for each m/z value.
% %The original reference for the alignment code is:
% Mantini,Petrucci,Pieragostino, Del Boccio, Di Nicola, Di Ilio, Federici, 
% Sacchetta, Comani and Urbani (2007). "LIMPIC: a computational method for 
% the separation of protein MALDI-TOF-MS signals from noise." 
% BMC Bioinformatics 8: 101.
[Peaks Intensity] = multiple_spectra_KL4(data,1,1);


%2. Find the elemental formulas for your peaks. We've provided our
%algorithm below (based on Kujawinski & Behn, Analytical Chemistry 2006)
%and modified in Kujawinski et al., GCA, 2009. 
%citation details: 
%Kujawinski, E. B. and M. D. Behn (2006). "Automated analysis of electrospray 
%ionization Fourier-transform ion cyclotron resonance mass spectra of natural 
%organic matter." Analytical Chemistry 78: 4363-4373.
%Kujawinski, E. B., K. Longnecker, N. V. Blough, R. Del Vecchio, L. Finlay, 
%J. B. Kitner and S. J. Giovannoni (2009). "Identification of possible source 
%markers in marine dissolved organic matter using ultrahigh resolution 
%mass spectrometry." Geochimica et Cosmochimica Acta 73: 4384-4399.

% Note that you will need to pick the right conversion to neutral masses 
%based on your ionization mode.


%2a. Convert to neutral mass (this assumes you calibrated on exact masses of
%ions, not on neutral masses).
H = 1.007825032;
elec = 5.4858e-4;
%IUPACpeak = Peaks + H - elec; % for negative ion mode
IUPACpeak = Peaks - H + elec; % for positive ion mode, Na-adducts are also likely

%2b. Call the formula determination algorithm
%you will need to change this next line to the location of the database
%(LongneckerKujawinski_fullCompoundList.mat)
fDir = 'C:\Documents and Settings\Krista\My Documents\MSdataAnalysis';
path(path,fDir);
load LongneckerKujawinski_fullCompoundList.2016.11.21.mat
[formulas elementOrder] = findformula_useList_KL17(IUPACpeak, zeros(size(IUPACpeak)), 1, 20, 500,fullCompoundList,'HAcap',1);

%3 Now go through and find the C13-containing formulas
FormulasC13 = quick13C_KL_1(formulas,IUPACpeak,0);
