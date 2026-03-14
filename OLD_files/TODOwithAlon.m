% locR2=locR/100;
% qrs_i_raw_ref2=qrs_i_raw_ref/2000;

% TODO order

final=zeros(length(qrs_i_raw_ref2),2);
for i=1:length(qrs_i_raw_ref2)
    diffI = abs(qrs_i_raw_ref2(i)-locR2);
    [val,loc]=min(diffI);
    final(i,1) = qrs_i_raw_ref2(i);
    final(i,2) = locR2(loc);
 
end

diffFinal = 60*diff(final)

diffOfTwo= diffFinal(:,1)-diffFinal(:,2);
meanOfTwo = mean(diffFinal,2)

figure;
scatter(meanOfTwo,diffOfTwo,'*')

figure; hold on;
plot(diffFinal(:,1),'.');
plot(diffFinal(:,2),'.');

finalMed=[medfilt1(diffFinal(:,1)) , medfilt1(diffFinal(:,2))]
figure; hold on;
plot(medfilt1(diffFinal(:,1),10));
plot(medfilt1(diffFinal(:,2),10));

diffOfMed= finalMed(:,1)-finalMed(:,2);
meanOfMed = mean(finalMed,2)

figure;
scatter(meanOfMed,diffOfMed,'*')

varofMed=cov(finalMed)
varNoMed=cov(diffFinal)

%find correlation between axes of finalMed(:1) and finalMed(:,2)
correlationMed = corr(finalMed(:,1), finalMed(:,2))
disp(['Correlation between the two axes of finalMed: ', num2str(correlationMed)]);

%same for diffFinal
% Calculate and display the correlation for diffFinal
correlationDiff = corr(diffFinal(:,1), diffFinal(:,2))
disp(['Correlation between the two axes of diffFinal: ', num2str(correlationDiff)]);
corrcoef



%%%% once this looks good, loop and save all data from all patients
% save into excel all data to see how our code is doing
% score ourselves
%  only then, consider kalman filter