fname = dir('Screen_data');
fname(1:2) = [];
fnum = numel(fname);
plate = cell(1,fnum);

    
%% Correlations between siRNA replicates

siRNA_corr = 1;
siRNA_null = 1;
Call = 0;
CallG = 0;
for i = 1:2:fnum
    i
    
    fname(i).name;
    
    D1 = readtable(['Screen_data/' fname(i).name]);
    D2 = readtable(['Screen_data/' fname(i+1).name]);
    
    mask = isnan(D1{:,6}) | isnan(D2{:,6});
    D1(mask,:) = [];
    D2(mask,:) = [];

    Dmat1 = D1{:,5:148};
    Dmat1 = zscore(Dmat1);
    
   
    Dmat2 = D2{:,5:148};
    Dmat2 = zscore(Dmat2);
    % ---- Look for plate effects ----
    %{
    subplot(2,4,(i+1)/2);
    
    [~,Dpc1] = pca(Dmat1);
    [~,Dpc2] = pca(Dmat2);
    
    P1 = [];
    P2 = [];
    
    for j = 1:max(D1{:,1})
        for k = 1:max(D1{:,2})
            v = find(j == D1{:,1} & k == D1{:,2});
            if v
                P1(j,k) = Dpc1(v,1);
                P2(j,k) = Dpc2(v,1);
            end
        end
    end
    imagesc([P1 P2])
    %}
    % ---- Test replicate correlation to screen correlation ----
    %{
    unique_siRNA = unique(D1{:,3});
    
    for j = 1:numel(unique_siRNA)
        mask = strcmp(D1{:,3},unique_siRNA(j));
       
        corr_genes = corr([Dmat1(mask,:); Dmat2(mask,:)]');
        siRNA_corr = [siRNA_corr corr_genes(:)'];
        
        R = randperm(size(Dmat1,1));
        corr_null1 = corr([Dmat1(mask,:); Dmat2(R(1:2),:)]');
        corr_null2 = corr([Dmat1(R(1:2),:); Dmat2(mask,:)]');
        siRNA_null = [siRNA_null corr_null1(:)' corr_null2(:)'];
    end
    
    %}
    %correlations between wells
    subplot(2,4,(i+1)/2);
    
    Mcorr = corr(Dmat1',Dmat2');
    Cp = Mcorr .* eye(size(Mcorr));
    Cn = Mcorr .* ~eye(size(Mcorr));
    
    Cp(Cp == 0) = [];
    Cn(Cn == 0) = [];
    
    Cpm = sum(Cp(:))/numel(Cp);
    Cnm = sum(Cn(:))/numel(Cn);
    
    Call = [Call Cp(:)' Cn(:)'];
    CallG = [CallG ones(1,numel(Cp(:)))*i ones(1,numel(Cn(:)))*(i+1)];
    
    image(corr(Dmat1',Dmat2')*64);
    
    title([fname(i).name ' ' num2str(Cpm) ' ' num2str(Cnm)]);
    
    %Mcorr = corr(Dmat1',Dmat2');
    
    % ---- Check wells match ----
    
    %for i = 1:min(size(Dmat1,1),size(Dmat2,1))
    %    a = a + strcmp(D1{i,3},D2{i,3});
    %end
    %size(Dmat2,1)
    
end
%%
figure
boxplot(Call,CallG)
    %%
siRNA_corr(siRNA_corr == 1) = [];
siRNA_null(siRNA_null == 1) = [];

hold on

histogram(siRNA_corr,0:0.05:1,'Normalization','pdf');
histogram(siRNA_null,0:0.05:1,'Normalization','pdf');