%%
while 1
    Coor = neuron.show_contours(0.6,[],[],1);
    inp=input('plot(p) / delete(d) / merge(m) / quit(q):','s');
    switch inp
        case 'p'
            ind=input('enter rois to PLOT:');
            figure; plot(neuron.C(ind,:)');
            pause;
        case 'd'
            ind=input('enter rois to DELETE:');
            neuron.delete(ind);
        case 'm'
            ind=input('enter two rois to MERGE:');
            neuron.manual_merge_multi_pairs({ind});
        case 'q'
            break;
    end
    delete(gcf);
end
