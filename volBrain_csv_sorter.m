clear
clc
male = importdata('bounds_male.csv'); 
female = importdata('bounds_female.csv');
filesToSort = dir('*report_p*'); %1x555
resFile = 'volBrain_Roger_results.txt'


for i = 1:length (filesToSort)
    fps= strsplit(filesToSort(i).name,'_');
    partID = fps{2};

    colCtr = 1;
    fnm = fullfile(filesToSort(i).folder, filesToSort(i).name);
    fnmData =importdata(fnm);
    sex = lower(fnmData.textdata{4});
    age = fnmData.textdata{6};
    fnmData.textdata = fnmData.textdata(:,9:end); %omit first 8 labels which are not needed
    for j = 1:(length(fnmData.data)) 
        %just gather percent (normed for size) and asymmetry scores!
        if  contains(fnmData.textdata{1,j},'%')  || contains(fnmData.textdata{1,j},'sym')
            matched(colCtr).name = fnmData.textdata{1,j};
            areas{colCtr} = fnmData.textdata{1,j};
            matched(colCtr).value = fnmData.data(j);
            colCtr = colCtr +1;
        end
    end

    %age gets columns for ROI {mean for age, min for age, max for age}
    if strcmp(sex,'male')
        selectedData = male.data(:,(str2num(age)*3+1:str2num(age)*3+3));
    else
        selectedData = female.data(:,(str2num(age)*3+1:str2num(age)*3+3));
    end
    
    %get some meaningful results
    format short g
    for k = 1:length(selectedData)

        perdat(i,k).partID = partID;
        perdat(i,k).area =  matched(k).name;
        perdat(i,k).subjScore = matched(k).value;
        perdat(i,k).low = selectedData(k,1);
        perdat(i,k).med = selectedData(k,2);
        perdat(i,k).high = selectedData(k,3);
        tempLine = [perdat(i,k).partID ' ' perdat(i,k).area  ' ' num2str(perdat(i,k).subjScore) ' ' num2str(perdat(i,k).low) ' ' num2str(perdat(i,k).med) ' ' num2str(perdat(i,k).high)];
    
        %measures of interest
        %less more or equal to mean
        if(perdat(i,k).subjScore < (perdat(i,k).med))
            perdat(i,k).moreless = '-1';
        elseif(perdat(i,k).subjScore > (perdat(i,k).med))
             perdat(i,k).moreless = '1';
        else
             perdat(i,k).moreless = '0';
             %fprintf  ('%d %d',i,k)
        end

        if(perdat(i,k).subjScore < (perdat(i,k).med - perdat(i,k).low))
               perdat(i,k).sigmoreless = '-1';
        elseif(perdat(i,k).subjScore > (perdat(i,k).med + perdat(i,k).high))
               perdat(i,k).sigmoreless = '1';
        else
                perdat(i,k).sigmoreless = '0';
        end
        
        n = perdat(i,k).subjScore - perdat(i,k).med;
        if n > 0
            d = (perdat(i,k).high - perdat(i,k).med);
        else
            d = (perdat(i,k).med - perdat(i,k).low);
        end
        perct= n/d;
        perdat(i,k).percMeas = perct * 100 ;
       
        %fprintf('%s\n',tempLine)

    end

end








%queryTerm = 'Postcentral gyrus total volume %'
% %pass queryTerm to match and get output of diff from zero and mean
% queryAreas(queryTerm,filesToSort,selectedData,perdat)

%ttests
for h = 1:length(areas)
    queryTerm = areas{1,h};
    queryAreasTtest(queryTerm,filesToSort,selectedData,perdat,resFile); 
end

%chisquares
for h = 1:length(areas)
    queryTerm = areas{1,h};
    chiSquare(queryTerm,filesToSort,selectedData,perdat,resFile); 
end

fclose('all')






function q = queryAreasTtest (queryTerm,filesToSort,selectedData,perdat,resFile) 
    clear returnList;
    for i = 1:length(filesToSort)
        for j = 1:length(selectedData)
            if (contains(perdat(i,j).area,queryTerm))
                returnList(i) =  perdat(i,j).percMeas;
            end
        end
    end
    returnList = returnList';
    format short g
    [h p] = ttest(returnList);
    m = mean(returnList);
    if (m > 0)
        direc = 'Greater in Athletes';
    elseif (m< 0) 
        direc = 'Lesser in Athletes';
    else
        direc = 'Void';
    end

    fprintf('PercMeas Diff for %s : p = %.3f, %s, mean = %.3f\n',queryTerm,round(p,3),direc,m);
    fid = fopen(resFile,'a');
    fprintf(fid,'PercMeas Diff for %s : p = %.3f, %s, mean = %.3f\n',queryTerm,round(p,3),direc,m);
    fclose(fid);

end

function q = chiSquare (queryTerm,filesToSort,selectedData,perdat,resFile)
    for i = 1:length(filesToSort)
        for j = 1:length(selectedData)
            if (contains(perdat(i,j).area,queryTerm))
                returnListmoreless(i) =  str2num(perdat(i,j).moreless);
                returnListsigmoreless(i) =  str2num(perdat(i,j).sigmoreless);   
            end 
        end
    end

    a = abs(sum(returnListmoreless(returnListmoreless<0)));
    b = abs(sum(returnListmoreless(returnListmoreless>0)));
    c = abs(sum(returnListsigmoreless(returnListsigmoreless<0)));
    d = abs(sum(returnListsigmoreless(returnListsigmoreless>0)));

    warning('off');
    contab = [27,28;a,b];
    p = chisquarecont(contab);
    if(a > b) 
        direc = 'Lesser in Athletes';
    elseif (b>a)
        direc = 'Greater in Athletes';
    else
        direc = "Void";
    end

    fprintf('ChiSqr Diff for %s : p = %.3f %s : Sig-%.0f Sig+%.0f\n',queryTerm,p,direc,c,d);
    fid = fopen(resFile,'a');
    fprintf(fid,'ChiSqr Diff for %s : p = %.3f %s : Sig-%.0f Sig+%.0f\n',queryTerm,p,direc,c,d);
    fclose(fid);

end
    



