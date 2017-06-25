function [A,C,C_raw,S] = prune_and_merge_neurons_postCNMFE(neuron)
figure;
subplot(10,1,1);
bar(std(neuron.C_raw')./std(neuron.C'));
%%
subplot(10,1,[2:10]);
plot_contours(neuron.A, reshape(sum(zscore(neuron.A(:,:)),2),240,320), 0.8, 1, [], [], 2);
%%
figure;
n_neurons = size(neuron.C,1);
locs = 1:n_neurons; % setdiff(1:n_neurons,[17 18 34 36 37 38]);
n_neurons = numel(locs);
for plotnum = n_neurons:-1:1
    h = axes;
    plot(neuron.C(plotnum,:));
    hold on;
    plot(neuron.C_raw(plotnum,:));
    h.XTick = [];
    h.YTick = [];
    pos = [0.1300    0.9213    0.7750    0.0037]; pos(2) = (plotnum-1)/n_neurons; pos(4) = 1/n_neurons;
    h.Position = pos;
    h.Box = 'off';
    h.Color = 'none';
    text(0,0,num2str(plotnum))
end
%% merging
flag = 1;
idxs = {};
while flag==1
    idx = input('Which neurons to merge? ([] to end):');
    if isempty(idx)
        flag = 0;
    else
        idxs = {idxs{:} idx};
    end
end
%% deleting
delidx = input('Which neurons to DELETE? ([] for none):');
%%
notidx = 1:n_neurons;
notidx = setdiff(notidx,unique([idxs{:}]));
notidx = setdiff(notidx,delidx);
A=[];C=[];C_raw=[];S=[];
for i = 1:numel(idxs)
    A = [A mean(neuron.A(:,idxs{i}),2)];
    C = [C; sum(neuron.C(idxs{i},:))];
    C_raw = [C_raw; sum(neuron.C_raw(idxs{i},:))];
    S = [S; sum(neuron.S(idxs{i},:))];
end
A = [A neuron.A(:,notidx)];
C = [C; neuron.C(notidx,:)];
C_raw = [C_raw; neuron.C_raw(notidx,:)];
S = [S; neuron.S(notidx,:)];

        