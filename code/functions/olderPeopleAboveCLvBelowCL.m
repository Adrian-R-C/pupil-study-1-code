function olderPeopleAboveCLvBelowCL()
% olderPeopleAboveCLvBelowCL is just a function that has been used to
% compare people above chance level vs. at chance level in the older group.
% It can be used to do the same in younger people as well.
%
%   Code created on May 7, 2025 by
%   Adrian RUIZ CHIAPELLO
%   Centre de Recherche Cerveau et Cognition
%   CNRS / Toulouse University.

cpt1=0;
oldRecognitionRates={app.recognitionRates{find(app.whosOld)}};

for cpt2=1:numel(app.medianPupilResponsesOld  )
    
    if ~ismember(cpt2, [10, 44])
        cpt1=cpt1+1;
        matrixAll(cpt1,:)=app.medianPupilResponsesOld{cpt2}.clean(:,5);
        recogAll(cpt1) = oldRecognitionRates{cpt2}(1);
        baselineAll(cpt1) = app.pupilMetricsOld(cpt2,5,1).meanBaselineValue;

    end
end
%Bien


%%
%%%ESSAYONS CHEZ LES JEUNES
% cpt1=0;
% for cpt2=1:numel(app.medianPupilResponsesYoung )
%     
%     if ~ismember(cpt2, [18, 23, 41, 44])
%         cpt1=cpt1+1;
%         matrixAllY(cpt1,:)=app.medianPupilResponsesYoung{cpt2}.clean(:,5);
%     end
% end
% %Bien
% 
% youngRecognitionRates={app.recognitionRates{find(app.whosYoung)}};
% cpt1=0;
% for cpt2=1:numel( youngRecognitionRates )
%     
%     if ~ismember(cpt2, [18, 23, 41, 44])
%         cpt1=cpt1+1;
%         recogAllY(cpt1) = youngRecognitionRates{cpt1}(1);
%     end
% end

%%%


matrixAboveCL=matrixAll(recogAll>=0.5900, :)';
matrixBelowCL=matrixAll(recogAll<0.5900, :)';

meanAboveCL=mean(matrixAboveCL,2);
meanBelowCL=mean(matrixBelowCL,2);

stdAboveCL=std(matrixAboveCL,0,2);
stdBelowCL=std(matrixBelowCL,0,2);

nAboveCL=size(matrixAboveCL,2);
nBelowCL=size(matrixBelowCL,2);
nSigma=app.sigmaForConfidenceInterval;


AboveCL_UPPER=meanAboveCL + nSigma * stdAboveCL / sqrt(nAboveCL);
AboveCL_LOWER=meanAboveCL - nSigma * stdAboveCL / sqrt(nAboveCL);

BelowCL_UPPER=meanBelowCL + nSigma * stdBelowCL / sqrt(nBelowCL);
BelowCL_LOWER=meanBelowCL - nSigma * stdBelowCL / sqrt(nBelowCL);

timepoints = (1:1500)/300;

%%

figure
plot(timepoints, matrixAboveCL, 'Color', app.colors(1,:), 'LineWidth', 1.2)
hold on
plot(timepoints, matrixBelowCL, 'Color', app.colors(2,:), 'LineWidth', 1.2)
legend( "Above CL (n=21)","At CL (n=25)" ) 

%%

figure
plot(timepoints,meanAboveCL, 'Color', app.colors(1,:) , 'LineWidth', 2)
hold on 
plot(timepoints,meanBelowCL, 'Color', app.colors(2,:), 'LineWidth', 2)
fill([timepoints, fliplr(timepoints)], [AboveCL_LOWER', fliplr(AboveCL_UPPER')], app.colors(1,:), 'FaceAlpha', 0.3, 'LineStyle', 'none');
fill([timepoints, fliplr(timepoints)], [BelowCL_LOWER', fliplr(BelowCL_UPPER')], app.colors(2,:), 'FaceAlpha', 0.3, 'LineStyle', 'none');
                        
[significance significanceCorrected clusterCharacteristics] = CBPT_ztValueThreshold(matrixAboveCL, matrixBelowCL, app.pThresholdCluster, app.numberOfPermutations, "independent", "both", "parametric");

plotGrayAreaSignificance(app, significance, -10, 20, timepoints)    

legend( "Above CL (n=21)","At CL (n=25)" ) 

[h p]=ttest2(matrixAboveCL', matrixBelowCL', "tail", "both")
figure
plot(timepoints, p)
%% Calcul métrique

for cpt=1:nAboveCL
    [peakDilationAboveCL(cpt) indiceMax(cpt)]=max(matrixAboveCL(:,cpt));
    latencyMax(cpt)=timepoints(indiceMax(cpt));
    meanDilationAboveCL(cpt)=mean(matrixAboveCL(:,cpt)>0);
    dilationVelocityAboveCL(cpt)=peakDilationAboveCL(cpt)/latencyMax(cpt);
end


for cpt=1:nBelowCL
    [peakDilationBelowCL(cpt) indiceMax(cpt)]=max(matrixBelowCL(:,cpt));
    latencyMax(cpt)=timepoints(indiceMax(cpt));
    meanDilationBelowCL(cpt)=mean(matrixBelowCL(:,cpt)>0);
    dilationVelocityBelowCL(cpt)=peakDilationBelowCL(cpt)/latencyMax(cpt);
end

meanBaselineValueAboveCL=baselineAll(recogAll>=0.5900)
meanBaselineValueBelowCL=baselineAll(recogAll<0.5900)

% t test indépendant 
% max dilatation
[h1 p1]=ttest2(peakDilationAboveCL,peakDilationBelowCL, "tail", "both");
% vitesse dilatation
[h2 p2]=ttest2(dilationVelocityAboveCL,dilationVelocityBelowCL, "tail", "both");
% dilatation moyenne
[h3 p3]=ttest2(meanDilationAboveCL,meanDilationBelowCL, "tail", "both");
% baseline moyenne
[h4 p4]=ttest2(meanBaselineValueAboveCL,meanBaselineValueBelowCL, "tail", "both");


%Plot de la baseline

donnees   = [meanBaselineValueAboveCL, meanBaselineValueBelowCL];
group_inx = [ones(1,numel(meanBaselineValueAboveCL)), 2.*ones(1,numel(meanBaselineValueBelowCL))];
condition_names=[{'Older above CL'},{'Older at CL'}];

c=[0.98,0.40,0.26 ; 0.25,0.00,0.66];

figure

hold on
h = daboxplot(donnees,'groups',group_inx,'mean',1,'color',c,...
    'xtlabels',condition_names, 'outliers', 0,...
        'scatter',2,'scattersize',100,'scatteralpha',0.6,'boxalpha',0.8);

ylabel('Performance');
set(gca,'FontSize',12)

% H=sigstar({[1,2]},[p4]);
axis square
%%


