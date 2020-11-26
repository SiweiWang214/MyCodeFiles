function [cfree,cknob,chromatri] = calchromknob(ring0,ind_sextfam)

for i = 1:length(ind_sextfam)
    pb{i} = atgetfieldvalues(ring0,ind_sextfam{i},'PolynomB');
    s0(i) = pb{i}{1}(3);
    for j = 1:10
        s_set(j,:) = s0(i) + 0.1*j;
        ring1 = atsetfieldvalues(ring0,ind_sextfam{i},'PolynomB',[0,0,s_set(j)]);
        [~,~,chrom(j,:)] = twissring(ring1,0,1:length(ring1)+1,'chrom');
    end
    coefx(i,:) = polyfit(s_set,chrom(:,1),1);
    coefy(i,:) = polyfit(s_set,chrom(:,2),1);
end

chromatri = [coefx(:,1),coefy(:,1)]';
cknob = pinv(chromatri);
[~,~,dmatri] = svd(chromatri);
cfree = dmatri(:,3:end);