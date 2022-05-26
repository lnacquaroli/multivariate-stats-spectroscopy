clc
clear
close all

%    ==========================================================
%    SOM_DEMO4 - DATA ANALYSIS USING THE SOM
%    ==========================================================

%    In this demo, the IRIS data set is analysed using SOM. First, the
%    data is read from ascii file (please make sure they can be found
%    from the current path), normalized, and a map is
%    trained. Since the data also had labels, the map is labelled.

sD = som_read_data('iris.data');     

sD = som_normalize(sD,'var');
sM = som_make(sD);
sM = som_autolabel(sM,sD,'vote');

%    VISUAL INSPECTION OF THE MAP
%    ============================

%    The first step in the analysis of the map is visual inspection.
%    Here is the U-matrix, component planes and labels (you may
%    need to enlarge the figure in order to make out the labels).

som_show(sM,'umat','all','comp',[1:4],'empty','Labels','norm','d');
som_show_add('label',sM.labels,'textsize',8,'textcolor','r','subplot',6);

%    From this first visualization, one can see that:
%     - there are essentially two clusters
%     - PetalL and PetalW are highly correlated
%     - SepalL is somewhat correlated to PetalL and PetalW
%     - one cluster corresponds to the Setosa species and exhibits
%       small petals and short but wide sepals
%     - the other cluster corresponds to Virginica and Versicolor
%       such that Versicolor has smaller leaves (both sepal and
%       petal) than Virginica
%     - inside both clusters, SepalL and SepalW are highly correlated

%    Next, the projection of the data set is investigated. A
%    principle component projection is made for the data, and applied
%    to the map. The colormap is done by spreading a colormap on the
%    projection. Distance matrix information is extracted from the
%    U-matrix, and it is modified by knowledge of zero-hits
%    (interpolative) units. Finally, three visualizations are shown:
%    the color code, with clustering information and the number of
%    hits in each unit, the projection and the labels.

echo off

%    From these figures one can see that: 
%     - the projection confirms the existence of two different clusters
%     - interpolative units seem to divide the Virginica
%       flowers into two classes, the difference being in the size of
%       sepal leaves

%    Finally, perhaps the most informative figure of all: simple
%    scatter plots and histograms of all variables. The species
%    information is coded as a fifth variable: 1 for Setosa, 2 for
%    Versicolor and 3 for Virginica. Original data points are in the
%    upper triangle, map prototype values on the lower triangle, and
%    histograms on the diagonal: black for the data set and red for
%    the map prototype values. The color coding of the data samples
%    has been copied from the map (from the BMU of each sample). Note
%    that the variable values have been denormalized.

echo off

%    This visualization shows quite a lot of information:
%    distributions of single and pairs of variables both in the data
%    and in the map. If the number of variables was even slightly
%    more, it would require a really big display to be convenient to
%    use.

%    From this visualization we can conform many of the earlier
%    conclusions, for example: 
%     - there are two clusters: 'Setosa' (blue, dark green) and 
%       'Virginica'/'Versicolor' (light green, yellow, reds)
%       (see almost any of the subplots)
%     - PetalL and PetalW have a high linear correlation (see
%       subplots 4,3 and 3,4)
%     - SepalL is correlated (at least in the bigger cluster) with
%       PetalL and PetalW (in subplots 1,3 1,4 3,1 and 4,1)
%     - SepalL and SepalW have a clear linear correlation, but it
%       is slightly different for the two main clusters (in
%       subplots 2,1 and 1,2)   

%    CLUSTERING OF THE MAP
%    =====================

%    Visual inspection already hinted that there are at least two
%    clusters in the data, and that the properties of the clusters are
%    different from each other (esp. relation of SepalL and
%    SepalW). For further investigation, the map needs to be
%    partitioned.

%    Here, the KMEANS_CLUSTERS function is used to find an initial
%    partitioning. The plot shows the Davies-Boulding clustering
%    index, which is minimized with best clustering.

subplot(1,3,1)
[c,p,err,ind] = kmeans_clusters(sM, 7); % find at most 7 clusters
plot(1:length(ind),ind,'x-')
[dummy,i] = min(ind)

cl = p{i};

%    The Davies-Boulding index seems to indicate that there are
%    two clusters on the map. Here is the clustering info
%    calculated previously and the partitioning result: 

subplot(1,3,2)
som_cplane(sM,Code,Dm)
subplot(1,3,3)
som_cplane(sM,cl)

%    You could use also function SOM_SELECT to manually make or modify
%    the partitioning.

%    After this, the analysis would proceed with summarization of the
%    results, and analysis of each cluster one at a time.
%    Unfortunately, you have to do that yourself. The SOM Toolbox does
%    not, yet, have functions for those purposes.

%    MODELING
%    ========

%    One can also build models on top of the SOM. Typically, these
%    models are simple local or nearest-neighbor models. 

%    Here, SOM is used for probability density estimation. Each map 
%    prototype is the center of a gaussian kernel, the parameters
%    of which are estimated from the data. The gaussian mixture
%    model is estimated with function SOM_ESTIMATE_GMM and the
%    probabilities can be calculated with SOM_PROBABILITY_GMM.

[K,P] = som_estimate_gmm(sM,sD);
[pd,Pdm,pmd] = som_probability_gmm(sD,sM,K,P);

%    Here is the probability density function value for the first data
%    sample (x=sD.data(:,1)) in terms of each map unit (m):

som_cplane(sM,Pdm(:,1))
colorbar
title('p(x|m)')

%    Here, SOM is used for classification. Although the SOM can be
%    used for classification as such, one has to remember that it does
%    not utilize class information at all, and thus its results are
%    inherently suboptimal. However, with small modifications, the
%    network can take the class into account. The function
%    SOM_SUPERVISED does this.

%    Learning vector quantization (LVQ) is an algorithm that is very
%    similar to the SOM in many aspects. However, it is specifically
%    designed for classification. In the SOM Toolbox, there are
%    functions LVQ1 and LVQ3 that implement two versions of this
%    algorithm.

%    Here, the function SOM_SUPERVISED is used to create a classifier
%    for IRIS data set:

sM = som_supervised(sD,'small');
som_show(sM,'umat','all');
som_show_add('label',sM.labels,'TextSize',8,'TextColor','r')

sD2 = som_label(sD,'clear','all'); 
sD2 = som_autolabel(sD2,sM);       % classification
ok = strcmp(sD2.labels,sD.labels); % errors
100*(1-sum(ok)/length(ok))         % error percentage (%)
