function [formula elementOrder] = findformula_useList_KL12(peak_mass, peak_int, formula_error, relation_error, mass_limit,fullCompoundList,sortType,showProgress)
%function [formula elementOrder] = findformula_useList_KL12(peak_mass, peak_int, formula_error, relation_error, mass_limit,fullCompoundList,sortType,showProgress)
%
%
%%%input the following variables:
% peak_mass - list of masses 
% peak_int - intensity of each peak (can currently be set to zeros)
%
% FORMULA_ERROR is the error for the elemental formula determination by the
% function, chemform, in ppm; 
%
% RELATION_ERROR is the window allowed % for identification of a relationship 
% between two peaks, in ppm.
%
% MASS_LIMIT is the maximum mass for the brute force calculation
%
% fullCompoundList - the database of possible compounds calculated by
% K Longnecker (LongneckerKujawinski_fullCompoundList.mat)
%
% sortType can be one of the following choices:
% 1) 'lowestSP' select the formula with the lowest # of S and P
% 2) 'HAcap' sorts based on the lowest number of S,P,N, and then only formulas
%       with P <= 1 or S <= 3 are considered valid
% 3) 'lowestError' to sort on the formula which has the lowest error from
%       the measured mass
% 
% This iteration % of the program uses only chemical relationships that 
% are common to refractory DOM such as humic acids. It does not include many
% biologically-relevant (i.e., metabolic) reactions.
% 
%%% The output is 
% formula : the elemental formulas for the peak list [C H O N C13 S P Na Error].
% elementOrder : a reminder of the order of elements by column
%
% original version of this algorithm published as: Kujawinski and Behn. 2006. 
% Automated analysis of electrospray ionization Fourier-transform ion 
% cyclotron resonance mass spectra of natural organic matter. 
% Analytical Chemistry 78:4363-4373. 
% Largest change has been to use a database to find the formulas rather
% than recalculating all of the possible formulas each time.
%
% Elizabeth Kujawinski Behn, May 2005
% Last updated: November 6, 2005
% updated version received by K Longnecker from LK 8/8/08
% KL changing relations to structures 9/12/08
%%KL added 9/12/08: put in a check to make sure that the format of the data 
%%are as needed for this function. peak_int and peak_mass both need to be 
%%multiple rows and one column
%%
%KL 1/2/09 - change the list sorting to consider both the minimum number of
%non-oxygen heteroatoms AND the lowest error - get around the case where
%can have different formulas with the same number of non-oxygen heteroatoms
%KL 1/6/09 - convert the results of getMATfile to double precision to get
%the right answer in the mass calculations! and added a line to Check7GR to
%require all elements to be positive (bc can't have negative number of
%elements in a compound, so why allow it)
%KL 1/17/09 findformula_uselist_KL7_SPsort - Still finding too much S and P
%- at the expense of N, try sorting to bias against S and P, but not N
%this involves changes both in the useMATfile function embedded and in the
%loops for building at the higher masses
%KL 1/20/09 - had a requirement that to keep a new peak, has to be within
%half of the given formula_error if there already is an old formula IF the
%formula is a comparison bw an existing formula found through the relations
%and a new formula found with the database (previously brute force)
%KL 4/15/09 - IF there is more than one formula, cap the number of S and P, 
%and then select based on lowest number of N, S, and P (sortType =
%'HA_cap') (will also test an overall cap on the number of S and P). Also,
%get rid of the relation which switches H for Na (will deal with that
%later)
%KL 4/30/09 - cleaning up, and changing the names of the sortTypes
%KL 5/1/09 - only change is to correct the neutral mass check bc sent out
%the wrong version yesterday
%KL 5/11/09 - changing the neutral mass check yet again
%KL 9/23/09 - change useMATfile to send out startform as double precision
%KL 9/1/2011 - correcting problem found by Dan Baluha (Univ. of Maryland)
%where in one of the calculations, the error calculation was multiplied by
%1e-6, rather than 1e6, therefore the error calculations were off by 1e-12.
%KL 4/22/2015 - Nikola Tolic (PNNL) found an error in the formula checking
%where we reference the wrong mass. v14 corrects that error

%check to be sure there are two variables going out
if nargout~=2
    error('Are you sure the right number of variables are going out?')
end


% check to make sure the masses are neutral masses
%KL 5/11/09 - change to if the masses have been corrected to neutral masses
%there should be more rounded numbers which are even numbers
ro = round(peak_mass); %first get the rounded masses
rm = mod(ro,2); %gets 1 = odd, 2 = even
even = length(find(rm == 0));
odd = length(find(rm == 1));
if odd > even
    error('Are you sure the masses have been converted to neutral masses?')
end
clear ro rm even odd

elementOrder = {'C' 'H' 'O' 'N' 'C13' 'S' 'P' 'Na' 'Error'};

relation_error = relation_error*1e-06;

% initial information about the elements
atomMasses %now this is an embedded function within findformula

%change 4/15/09 to not allow the Na/H switch
relations = [1 2 0 0 0 0 0 0; 1 4 -1 0 0 0 0 0; 0 2 0 0 0 0 0 0; ...
         2 4 1 0 0 0 0 0; 1 0 2 0 0 0 0 0; 2 2 1 0 0 0 0 0; 0 0 1 0 0 0 0 0;...
         1 1 0 0 0 0 0 0; 0 1 0 1 0 0 0 0;  ...
         0 0 3 0 0 0 1 0]; % matrix that defines each mass change
elements = [C; H; O; N; C13; S; P; Na]; % matrix that defines the elements of interest
nElements = length(elements);
grp = relations*elements;

% Find relations between peaks
f = length(peak_mass);

%KL: setup a structure rather than the multiple Grpdiff# matrices
GrpdiffK = struct('mass',num2cell(peak_mass),'grpdiff',[]);

for p=1:f
    for x=1:length(grp)
        %KL note: look for differences, and keep those less than the defined relation_error
        [temp] = find(abs((peak_mass(p)-peak_mass)/grp(x) - round((peak_mass(p)-peak_mass)/grp(x)))<=relation_error ...
            & peak_mass(p)~=peak_mass); % asks if any relations exist for functional group x
       
        if isempty(temp)==0  % some relations have been found
            %9/12/08 KL changing to use the structure
            GrpdiffK(p).grpdiff{x} = [peak_mass(p) peak_int(p); peak_mass(temp) peak_int(temp)];            
        end % ends if isempty loop
    end % ends for x loop
end % ends for p loop
% fprintf('Finished relations\n')

% Find formulas for peaks
formula = zeros(f,length(elements)+1);
for p=1:f   
    %y = isempty(GrpdiffK(p).grpdiff); KL: change this bc I keep forgetting what y is...
    if peak_mass(p) < mass_limit % checks mass, if < mass_limit      
        
        if isempty(GrpdiffK(p).grpdiff) % Grpdiff is empty but mass < mass_limit; KL note - so do the brute force formula finding           
            startform = useMATfile_KL4(peak_mass(p),fullCompoundList,formula_error);
 
            if startform(1:nElements) >= 0 %KL note: i.e. have to have positive numbers for at least some of the elements             
                tempmass = startform*elements; % find the theoretical mass
                Terror = abs((peak_mass(p)-tempmass)/tempmass)*1e06; % find the error: KL change to Terror bc 'error' is a matlab command
                formula(p,:) = [startform Terror]; % insert identified formula into master matrix
            end
            
        elseif ~isempty(GrpdiffK(p).grpdiff) % do this if Grpdiff is not empty (have relations)
            if formula(p,1)==0 % if formula is not known...
                startform = useMATfile_KL4(peak_mass(p),fullCompoundList,formula_error);
                               
                if startform(1:nElements) >= 0 
                    tempmass = startform*elements; 
                    Terror = abs((peak_mass(p)-tempmass)/tempmass)*1e06; 
                    formula(p,:) = [startform Terror]; 
                end
            elseif formula(p,1) > 0 % if formula is known already...
                %first get the brute force formula from the table
                tempform = useMATfile_KL4(peak_mass(p),fullCompoundList,formula_error);
                tempmass = tempform*elements;
                errornew = abs((peak_mass(p)-tempmass)/tempmass*1e06);
                
                check = Check7GR_KL2(formula(p,1:nElements),peak_mass(p)); %this seems redundant...but leave in
                
                %then decide whether to keep the brute force formula or the
                %formula which is already present (from the relations)
                if check == 0
                    formula(p,:) = [tempform errornew];
                elseif tempform(1:7) >= 0 & ...
                        (tempform(6)+tempform(7)) < (formula(p,6)+formula(p,7)) & ...
                        errornew < formula_error/2; %KL adding this last line 1/20/09 - only keep brute force if 
                                                    %it 'really' improves the formula (much smaller error)
                        
                    formula(p,:) = [tempform errornew];
                end
            end % ends formula loop
            
            %now go through Grpdiff and assign those formulas...so using
            %the formula decided to be ok - build from that using the relations
            startform = formula(p,1:length(elements));

            %%KL change to only do this for the places where Grpdiff is not empty 
            ge = cellfun(@isempty,GrpdiffK(p).grpdiff);
            gek = find(ge == 0);
           
            for x = gek
                ftemp = GrpdiffK(p).grpdiff{x}; %KL changing to the array
                [r,c] = size(ftemp);
                for t=2:r % iterate through related peaks, KL note: start with t = 2 bc t = 1 is the test mass 
                    numGps = round((ftemp(t,1)-ftemp(1,1))/grp(x,:)); % find # of grps in relation
                    if numGps > 0 % checks to be sure there is at least one GRP difference between two peaks
                        relpk = find(ftemp(t,1)==peak_mass);  % find related peak in original peak list
                        newform = startform + numGps*relations(x,:); % find formula for related peak
                        newmass = newform*elements;
                        errornew = abs((peak_mass(relpk)-newmass)/newmass)*1e06; % calculate error

                        checknew = Check7GR_KL2(newform,peak_mass(relpk)); 
                        checkold = Check7GR_KL2(formula(relpk,1:nElements),peak_mass(relpk));
                        
                        if checknew==1
                            if checkold==1 & ...
                                    (newform(6)+newform(7)) < (formula(relpk,6)+formula(relpk,7)) &...
                                    errornew <= formula_error
                                formula(relpk,:) = [newform errornew];
                            elseif checkold == 0  & errornew <= formula_error 
                                formula(relpk,:) = [newform errornew];
                            end % end formula check loop
                        end % ends positive formula loop
                        clear checknew checkold
                        
                    end % ends numGps loop
                end % ends t loop "rows within each depth"
            end % ends x loop "depths"
        end % ends y loop "non-empty Grpdiff"

       
    else % if MW > 500...KL note: (MW is the mass_limit set as input argument)
        if ~isempty(GrpdiffK(p).grpdiff) % is Grpdiff empty?  
            if formula(p,1)> 0 % if the formula is already known...
                startform = formula(p,1:nElements);
                %%KL change to only do this for the places where Grpdiff is not empty 
                ge = cellfun(@isempty,GrpdiffK(p).grpdiff);
                gek = find(ge == 0);

                for x = gek             
                    ftemp = GrpdiffK(p).grpdiff{x}; %KL changing to the array
                    [r,c] = size(ftemp);
                    
                    for t=2:r
                        numGps = round((ftemp(t,1)-ftemp(1,1))/grp(x,:));                        
                        relpk = find(ftemp(t,1)==peak_mass);
                        [checkold gf] = Check7GR_KL2(formula(relpk,1:nElements),peak_mass(relpk));
                        newform = startform + numGps*relations(x,:);
                        [checknew gf] = Check7GR_KL2(newform,peak_mass(relpk));
                        errornew = abs((peak_mass(relpk)-(newform*elements))/(newform*elements))*1e06;
                          
                        if checknew==1 & abs(numGps) <= 5 & errornew <= formula_error
                            if checkold==1 & ...
                                      (newform(6)+newform(7)) < (formula(relpk,6)+formula(relpk,7))                        
                                formula(relpk,:) = [newform errornew];
                            elseif checkold==0
                                formula(relpk,:) = [newform errornew];
                            end % end formula check loop
                        end % ends positive formula loop
                        clear checkold checknew
    
                    end % ends t loop
                end % ends x loop
                    
            else % if formula is not known - check for formulas in lower masses
                ge = cellfun(@isempty,GrpdiffK(p).grpdiff);
                gek = find(ge == 0);
                for x = gek 
                    ftemp = GrpdiffK(p).grpdiff{x}; %KL changing to the array
                    [r,c] = size(ftemp);
                    
                    for t=2:r
                        if ftemp(t,1) < ftemp(1,1)
                            low_m = find(ftemp(t,1)==peak_mass);
                            if formula(low_m,1) > 0
                                startform = formula(low_m,1:length(elements));
                                numGps = round(abs(ftemp(t,1)-ftemp(1,1))/grp(x,:));
                                newform = startform + numGps*relations(x,:);                                  
                                checknew = Check7GR_KL2(newform, peak_mass(p)); %%KL corrected 4/22/2015      
                                %KL corrected the following line 8/31/2011...note the sign of 1e6
                                %errornew = abs((peak_mass(p)-newform*elements)/(newform*elements))*1e-06; 
                                errornew = abs((peak_mass(p)-newform*elements)/(newform*elements))*1e06;
                                if abs(numGps)<= 5 & checknew==1 & errornew <= formula_error 
                                    formula(p,:) = [newform errornew];
                                end % ends formula check loop
                                clear checknew                              
                            end % ends finding lower formula
                        end % ends ftemp
                    end % ends t loop
                end % ends x loop
            end % ends formula loop
        end % ends GrpdiffK loop (was called the y loop, but I (KL) keep having to look up 'y')
    end % ends high mass loop
    if showProgress && ~rem(p,250)
        fprintf('Finished Peak %1.0f of %1.0f\n',p,f)
    end
    
end % ends p loop

% Kendrick mass stuff
KenMass = peak_mass/(C+2*H)*14;
NomIUPACMass = floor(peak_mass);
KMD = floor((NomIUPACMass-KenMass)*1000);
Zstar = mod(NomIUPACMass, 14)-14;
Kendrick_matrix = [peak_mass KenMass KMD Zstar];

% Kendrick loop
for p=1:f
    if formula(p,1) > 0
        startform = formula(p,1:nElements);
        Ktemp = Kendrick_matrix(find((Kendrick_matrix(:,3) == Kendrick_matrix(p,3)) & (Kendrick_matrix(:,4) ...
            == Kendrick_matrix(p,4))),:);
        [r,c] = size(Ktemp);
            
        for t=2:r
            numCH2 = floor((Ktemp(t,1)-Ktemp(1,1))/14);
            relpk = find(Ktemp(t,1)==peak_mass);
            CH2 = zeros(1,nElements);
            CH2(1:2) = [1 2]; %KL change this to be more general
            
            newform = startform + numCH2*CH2;
            checknew = Check7GR_KL2(newform,peak_mass(relpk));
            checkold = Check7GR_KL2(formula(relpk, 1:nElements),peak_mass(relpk));
            errornew = abs((peak_mass(relpk)-(newform*elements))/(newform*elements))*1e06;
            
            if checkold==1 & checknew==1 & errornew <= formula_error/2 & ...
                    (newform(6)+newform(7)) < (formula(relpk,6)+formula(relpk,7))
                formula(relpk,:) = [newform errornew];
            elseif checkold == 0 & checknew==1 & errornew <= formula_error/2
                formula(relpk,:) = [newform errornew];
            end % end formula check loop
            clear checknew checkold
         
        end % ends t loop
    end % ends formula loop
end % ends p loop  


%%%start the embedded functions here
%%%start the embedded functions here
    function atomMasses
    % initial information
    H = 1.007825032;
    C = 12;
    O = 15.994914622;
    N = 14.003074005;
    Na = 22.989767;
    Br = 79.904;
    Cl = 35.4527;
    C13 = 13.003354826;
    N15 = 15.000108;
    P = 30.973762;
    S = 31.972070;
    end %end atomMasses embedded function


    function [good goodformulas] = Check7GR_KL2(formulas, mass)
    % [GOOD, GOODFORMULAS] = Check7GR(FORMULAS, MASS) where GOOD is 1 if formula adheres to rules and
    % 0 if not; FORMULAS is a list of formulas that have been proposed by the
    % relationship algorithm (usually just one) and MASS is the neutral mass
    % (corrected from observed)

    % this function should check any formulas assigned by the relationships
    % with the 7-golden rules (BMC Bioinformatics 2007 8:105)

    % written August 2008 LK
    %KL modify 9/18/08 to go back to presenting good and goodformula
    %KL 1/8/09 - adding a check to make sure all the elements are positive
    %avoid having negative number of elements

    LoadAtomMasses 

    elements = [C; H; O; N; C13; S; P; Na];
    valence = [4; 1; 2; 3; 4; 2; 3; 1];

    keep = formulas;
    None = [-100, zeros(1,length(elements))];
    lowlimit = [1,1, zeros(1,length(elements)-2)];

    if keep(1)<lowlimit(1) & keep(2)<lowlimit(2)
        good = 0;
        goodformulas = None;
        return
    end

    % step 1: check number of elements possible within mass range
    if mass<=500
        uplimit = [39, 72, 20, 20, 0, 10, 9, 0];
        keep0 = keep(find((keep(:,1)+keep(:,5))<uplimit(1) & ...
            keep(:,2)<uplimit(2) & keep(:,3)<uplimit(3) & ...
            keep(:,4)<uplimit(4) & keep(:,6)<uplimit(6) & ...
            keep(:,7)<uplimit(7)),:);
    elseif mass>500 & mass<=1000
        uplimit = [78, 126, 27, 25, 0, 14, 9, 0];
        keep0 = keep(find((keep(:,1)+keep(:,5))<uplimit(1) & ...
            keep(:,2)<uplimit(2) & keep(:,3)<uplimit(3) & ...
            keep(:,4)<uplimit(4) & keep(:,6)<uplimit(6) & ...
            keep(:,7)<uplimit(7)),:);
    elseif mass>1000 & mass<=2000
        uplimit = [156, 180, 63, 32, 0, 14, 9, 0];
        keep0 = keep(find((keep(:,1)+keep(:,5))<uplimit(1) & ...
            keep(:,2)<uplimit(2) & keep(:,3)<uplimit(3) & ...
            keep(:,4)<uplimit(4) & keep(:,6)<uplimit(6) & ...
            keep(:,7)<uplimit(7)),:);
    end

    if isempty(keep0)==1
        good = 0;
        goodformulas = None;
        return
    end

    % step 2: check for valence rules
    oddVal = [keep0(:,2), keep0(:,4)];
    keep1 = keep0(find(rem(sum(oddVal,2),2)==0),:);

    if isempty(keep1)==1
        good = 0;
        goodformulas = None;
        return
    end

    % sum of valences greater than or equal to twice the maximum valence
    sumVal = keep1(:, 1:length(elements))*valence;
    uniVal = keep1(:,1:length(elements));
    uniVal(uniVal~=0)=1;
    [r,c] = size(keep1);
    for x=1:r
        maxVal(x) = max(uniVal(x, 1:length(elements)).*valence');
    end
    keep2 = keep1(find(sumVal>=2*maxVal'),:);
    clear sumVal

    if isempty(keep2)==1
        good = 0;
        goodformulas = None;
        return
    end

    % sum of valences greater than or equal to 2*atom#-1
    sumVal = keep2(:,1:length(elements))*valence;
    test = 2*sum(keep2(:,1:length(elements)),2)-1;
    keep3 = keep2(find(sumVal>=test),:);

    if isempty(keep3)==1
        good = 0;
        goodformulas = None;
        return
    end

    % step 3: check elemental ratios

    HC = keep3(:,2)./(keep3(:,1)+keep3(:,5));
    NC = keep3(:,4)./(keep3(:,1)+keep3(:,5));
    OC = keep3(:,3)./(keep3(:,1)+keep3(:,5));
    PC = keep3(:,7)./(keep3(:,1)+keep3(:,5));
    SC = keep3(:,6)./(keep3(:,1)+keep3(:,5));
    keep4 = keep3(find(HC>=0.2 & HC <= 3.1 & NC>=0 & NC<=1.3 & OC>=0 & OC<=1.2 &...
        PC>=0 & PC<=0.3 & SC>=0 & SC<=0.8),:);

    if isempty(keep4)==1
        good = 0;
        goodformulas = None;
        return
    end

    % step 4: check heteroatom abundances
    checkNOPS = find(keep4(:,3)>=1 & keep4(:,4)>=1 & keep4(:,6)>=1 & keep4(:,7)>=1);
    checkNOP = find(keep4(:,3)>=3 & keep4(:,4)>=3 & keep4(:,6)==0 & keep4(:,7)>=3);
    checkOPS = find(keep4(:,3)>=1 & keep4(:,4)==0 & keep4(:,6)>=1 & keep4(:,7)>=1);
    checkNPS = find(keep4(:,3)==0 & keep4(:,4)>=1 & keep4(:,6)>=1 & keep4(:,7)>=1);
    checkNOS = find(keep4(:,3)>=6 & keep4(:,4)>=6 & keep4(:,6)>=6 & keep4(:,7)==0);
    [row,col]=size(keep4);
    allr = [1:1:row]';
    check = sort([checkNOPS; checkNOP; checkOPS; checkNPS; checkNOS]);
    kp = [];
    for c=1:length(allr)
        tp = find(check==allr(c));
        if isempty(tp)==1
            kp = [kp; c];
        end
    end
    if isempty(checkNOPS)==0
        kp1 = find(keep4(checkNOPS,3)<20 & keep4(checkNOPS,4)<10 & keep4(checkNOPS,6)<3 & keep4(checkNOPS,7)<4);
        kp = [kp; checkNOPS(kp1)];
    elseif isempty(checkNOP)==0
        kp2 = find(keep4(checkNOP,3)<22 & keep4(checkNOP,4)<11 & keep4(checkNOP,7)<6);
        kp = [kp; checkNOP(kp2)];
    elseif isempty(checkOPS)==0
        kp3 = find(keep4(checkOPS,3)<14 & keep4(checkOPS,6)<3 & keep4(checkOPS,7)<3);
        kp = [kp; checkOPS(kp3)];
    elseif isempty(checkNPS)==0
        kp4 = find(keep4(checkNPS,4)<4 & keep4(checkNPS,6)<3 & keep4(checkNPS,7)<3);
        kp = [kp; checkNPS(kp4)];
    elseif isempty(checkNOS)==0
        kp5 = find(keep4(checkNOS,3)<14 & keep4(checkNOS,4)<19 & keep4(checkNOS,6)<8);
        kp = [kp; checkNOS(kp5)];
    end
    
    kp = sort(kp);
    keep5 = keep4(kp,:);

    if isempty(keep5)==1
        good = 0;
        goodformulas = None;
    return
    end

    %KL 1/8/09 - checking to make sure everything is a positive number
    %(i.e. can't have a negative number of elements in a formula)
    if keep5>=0;
        good = 1;
        goodformulas = keep5;
    else 
        good = 0;
        goodformulas = None;
    end
    

    %%LoadAtomMasses as a function within Check7GR
        function LoadAtomMasses
        H = 1.007825032;
        C = 12;
        O = 15.994914622;
        N = 14.003074005;
        S = 31.972070;
        P = 30.973762;
        Cl = 35.4527;
        Br = 79.904;
        Na =22.989767;
        C13 = 13.003354826;
        end

    end %end Check7GR as a function


    function elemformula = useMATfile_KL4(peak_mass, fullCompoundList, sWin)
    %function elemformula = useMATfile_KL4(peak_mass, fullCompoundList, sWin)
    %using the master list of compounds - can keep them in the array and then
    %use arrayfun to search through the different cells in sData
    %KL 12/15/08
    %1/2/09 - can have the case where multiple compounds have the same number
    %of non-oxygen heteroatoms, so that isn't always the best way to sort. If
    %that happens, from the options with the lowest number of non-oxygen
    %heteroatoms, chose the one with the lowest error
    %KL 1/6/09 - convert elemformula to double to get the precision I need
    %KL 9/23/09 - make sure to send out elemformula as double precision


    %decide how broad the m/z range should be - using sWin
    upperMZ = peak_mass + (peak_mass)*(sWin/1e6);
    lowerMZ = peak_mass - (peak_mass)*(sWin/1e6);

    findMass = arrayfun(@(x) find(x.mass > lowerMZ & x.mass < upperMZ), fullCompoundList,'uniformOutput',false);

    findMassIdx = cellfun(@isempty,findMass); %will spit out 1 if empty - actually want the 0 here
    kB = find(findMassIdx==0); %figure out where we need to go into sData - can be more than one cell

    testForm = [];
    testMass = [];
    testError = [];

    for a = 1:length(kB);
        %for easier use - pull the index into sData first
        idx = findMass{kB(a)};

        %then use that to get the testMass and testForm
        testMass = [testMass ; fullCompoundList(kB(a)).mass(idx)];
        testForm = [testForm ; fullCompoundList(kB(a)).formula(idx,:)];
        clear idx
        %in order to sort based on error (testing), need to actually
        %calculate the error
        testError = [testError ; abs((peak_mass - testMass)./testMass)*1e6];
    end
    clear a   
        
    if isempty(testForm)
        n = size(fullCompoundList(1).formula,2);
        elemformula = zeros(1,n); clear n
        elemformula(1) = -100;
    else       
        %%need to select from the final list:
        
        switch sortType
            case 'lowestSP' 
                %select the first elemental formula on the list (lowest # of
                %non oxygen heteroatoms

                %NonOxyHeteroAtoms = testForm(:,4)+testForm(:,6)+testForm(:,7);
                %1/20/09: change to test if only sort based on S and P
                NonOxyHeteroAtoms = testForm(:,6) + testForm(:,7);
                [iy,iz]=sort(NonOxyHeteroAtoms);

                %KL 1/2/09 addition - if multiple with the same, low, # of
                %non-oxygen heteroatoms, sort based on the lowest error AND
                %lowest number of non-oyxgen heteroatoms (S and P only here).
                k = find(NonOxyHeteroAtoms == iy(1)); %iy(1) will be the lowest # of non-oxy heteroatoms
                if length(k)>1;
                    %only consider the formulas with the low # non-oxy HA
                    smallForm = testForm(k,:);
                    smallError = testError(k,1);

                    %sort those based on the error:
                    [Y I] = sort(smallError);
                    elemformula = smallForm(I(1),:);

                else
                    %easy - can use the first on the list
                    elemformula = testForm(iz(1),:);    
                end
            
            case 'lowestError' 
                %sort based on the lowest error - for comparison
                [Y I] = sort(testError);
                elemformula = testForm(I(1),:); 
            
            case 'HAcap'
                %added 4/15/09 by KL cap the number of S and P atoms, after
                %selecting based on the lowest number of N, S, and P
                
                %calculate the # of non-oxygen heteroatoms
                NonOxyHeteroAtoms = testForm(:,4)+testForm(:,6)+testForm(:,7);
                [iy,iz]=sort(NonOxyHeteroAtoms);

                %KL 1/2/09 addition - if multiple with the same, low, # of
                %non-oxygen heteroatoms, need to do more sorting
                k = find(NonOxyHeteroAtoms == iy(1)); %iy(1) will be the lowest # of non-oxy heteroatoms
                if length(k)>1;
                    %only consider the formulas with the low # non-oxy HA
                    smallForm = testForm(k,:);
                    smallError = testError(k,1);

                    %then only consider formulas with P <= 1 or S <= 3
                    kPS = find(smallForm(:,7) <= 1 & smallForm(:,6) <= 3);                  
                    if length(kPS) > 1;
                        %if more than one formula passes the PS cap, sort
                        %those that do based on the smallest error
                        evenSmallerForm = smallForm(kPS,:);
                        evenSmallerError = smallError(kPS,:);
                        [Y I]= sort(evenSmallerError);
                        elemformula = evenSmallerForm(I(1),:);
                    elseif length(kPS) == 1
                        %in this case, there is only one formula left - use
                        %that one
                        elemformula = smallForm(kPS,:);
                    elseif isempty(kPS)
                        %nothing on this list passes the PS cap, so go back
                        %to the bigger list, and get the formula with the
                        %lowest error and the lowest # of N,S, and P
                        
                        %sort those based on the error:
                        [Y I] = sort(smallError);
                        elemformula = smallForm(I(1),:);
                    end
                    clear smallForm smallError evenSmallerForm evenSmallerError kPS                   
                else
                    %easy - can use the first on the list
                    elemformula = testForm(iz(1),:);    
                end
            otherwise %added this 9/23/09 - make it easier in case of typos
                error('sortType is not recognized...something is wrong')

        end
    end
    %startform comes out as single to save space, convert to double
    %to get the correct masses when doing calculations
    elemformula = double(elemformula);


    end %end useMATfile_KL4 as function

    
end %end findformula function