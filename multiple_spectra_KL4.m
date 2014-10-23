function [Peaks Intensity]=multiple_spectra_KL4(dataIn,ppm,showProgress)
%function [Peaks Intensity]=multiple_spectra_KL4(dataIn,ppm,showProgress)
%
% input:
% dataIn is the unaligned data - in a cell where each cell has 2 columns
% which are [peaks intensity]
% ppm is the maximum width for peak clustering (in parts per million)
% showProgress = 1 indicates you want to see where the algorithm is (more or
% less); 0 indicates the code should not show its iterations
% 
% output:
% Peaks is the vector containing the m/z of the peak classes
% Intensity is the matrix of peak amplitudes for the analyzed spectra
% 
%KL 11/21/08 using the LIMPIC project's alignment tool to align the different
%spectra - don't need the peak picking part of that algorithm. Really the
%only changes have been how the data go in and come out, though I did add
%comments in here bc there were none in the original m-file
%
%1/15/09 - version KL4 just changes to send out Peaks and Intensity
%with the peaks in one column - so that each sample is a different column
%in Intensity and Peaks is a [r x 1] matrix
%
% %The original reference is:
% Mantini,Petrucci,Pieragostino, Del Boccio, Di Nicola, Di Ilio, Federici, 
% Sacchetta, Comani and Urbani (2007). "LIMPIC: a computational method for 
% the separation of protein MALDI-TOF-MS signals from noise." 
% BMC Bioinformatics 8: 101.
%
%%contact information:
% Krista Longnecker
% Department of Marine Chemistry and Geochemistry
% Woods Hole Oceanographic Institution
% klongnecker@whoi.edu

data =[];
for a=1:length(dataIn)
    peak_mass = dataIn{a}(:,1);
    peak_intensity = dataIn{a}(:,2);
    ntp=length(peak_intensity);
    data=[data ; [peak_mass' ; peak_intensity'; a*ones(1,ntp)]'];
    clear ntp peak_intensity peak_mass k td
end
clear a

nspectra = length(unique(data(:,3)));

% creation of data matrices for peak clustering
[dataord(:,1),order]=sort(data(:,1));
dataord(:,2:3)=data(order,2:3);

Peaks=dataord(:,1)';
npoints=length(Peaks);
lista_nel=ones(1,npoints);
Intensity=zeros(nspectra,npoints);
for i=1:npoints
    pat=dataord(i,3);
    amp=dataord(i,2);
    Intensity(pat,i)=amp;
end


% iterative peak clustering
i=2;

while i < length(Peaks)
    if rem(i,10)==1 & showProgress %KL added this loop
        fprintf('iteration %d for %d peaks\n',i,length(Peaks)-1) 
    end
    
    value_b=Peaks(i-1);
    value=Peaks(i);
    value_p=Peaks(i+1);
    db=value-value_b;
    dp=value_p-value;
    thres=value*ppm/1000000;
    pat=find(Intensity(:,i) > 0); 
    logicval=true;
    
    %KL note: looking at the peak on the lower end
    for k=1:length(pat)
        logic =(Intensity(pat(k), i-1 ) == 0) ;
        logicval=logicval & logic;
    end

    if logicval == false
        db = inf;
    end
    logicval=true;
    
    %KL note: now do the same for the peak on the upper end
    for k=1:length(pat)
        logic =(Intensity(pat(k), i+1 ) == 0) ;
        logicval=logicval & logic;
    end
    
    if logicval == false
        dp = inf;
    end
    
    if db <= dp & db < thres    
      %KL note: peak1 and peak2 are closer together.
      %KL note: Take the average m/z of peak1 and peak2 and put that
      %into Peaks
      Peaks(i-1)=(lista_nel(i)*value+lista_nel(i-1)*value_b)/(lista_nel(i-1)+lista_nel(i));
      lista_nel(i-1)=lista_nel(i-1)+lista_nel(i); 
      %KL note: adding the intensity of peak1 and peak2
      Intensity(:,i-1)=Intensity(:,i-1)+Intensity(:,i);
      %and delete the index
      Peaks(i)=[];
      lista_nel(i)=[];
      Intensity(:,i)=[];
    elseif db > dp & dp < thres 
      %%KL note: peak2 and peak3 are closer together (bc distance between peak1 and
      %peak2 is greater)
      %KL note: Take the average m/z of peak2 and peak3 and put that
      %into Peaks
      Peaks(i+1)=(lista_nel(i)*value+lista_nel(i+1)*value_p)/(lista_nel(i+1)+lista_nel(i));
      lista_nel(i+1)=lista_nel(i+1)+lista_nel(i); 
      %KL note: adding the intensity of peak2 and peak3
      Intensity(:,i+1)=Intensity(:,i+1)+Intensity(:,i);
      %KL note: delete the index...
      Peaks(i)=[];
      lista_nel(i)=[];
      Intensity(:,i)=[];  
    else
        %KL note: otherwise, the distance bw the two peaks is greater than 
        %the threshold - move on to the next set
      i=i+1;
    end
end


lista_min=Peaks*(1-ppm/1000000); 
lista_max=Peaks*(1+ppm/1000000);

%how many peaks did we come up with?
nclassi=length(Peaks);

for i=2:nclassi
    if lista_min(i) < lista_max(i-1)
        lista_min(i)=mean(Peaks(i-1:i));
        lista_max(i-1)=lista_min(i);
    end
end

%KL 1/15/09 - change for version KL4
Peaks = Peaks';
Intensity = Intensity';