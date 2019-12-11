function  psth_out=psth_freq_fun(i,sponts,tabss,trels,nrep,dat,psth_freq,reptime, ts,Cohc,Cihc,binw,len,binwidth,B,A)

    sponts_concat = [sponts.LS(i,1:nrep(1)) sponts.MS(i,1:nrep(2)) sponts.HS(i,1:nrep(3))];
    tabss_concat = [tabss.LS(i,1:nrep(1)) tabss.MS(i,1:nrep(2)) tabss.HS(i,1:nrep(3))];
    trels_concat = [trels.LS(i,1:nrep(1)) trels.MS(i,1:nrep(2)) trels.HS(i,1:nrep(3))];


    %vihc = model_IHC_BEZ2018(pin,CF,nrep,dt,2*T,cohc,cihc,species);
    vihc_temp = model_IHC_BEZ2018(dat,psth_freq(i),1,ts,reptime,Cohc(i),Cihc(i),2);
    
    for ma = 1:nrep(1)
       spont = sponts_concat(ma);
       tabs = tabss_concat(ma);
       trel = trels_concat(ma);
       [psth500k_a_temp,~,~,~,~,~] = model_Synapse_BEZ2018(vihc_temp,psth_freq(i),1,ts,1,0,spont,tabs,trel);
        psth500k_a_single_fibers(ma,:) = psth500k_a_temp;
    end

    psth500k_a = sum(psth500k_a_single_fibers);
    
    for mb = 1:nrep(2)
       spont = sponts_concat(mb+nrep(1));
       tabs = tabss_concat(mb+nrep(1));
       trel = trels_concat(mb+nrep(1));
       [psth500k_b_temp,~,~,~, ~,~] = model_Synapse_BEZ2018(vihc_temp,psth_freq(i),1,ts,1,0,spont,tabs,trel);
       psth500k_b_single_fibers(mb,:) = psth500k_b_temp;
    end
 
    psth500k_b = sum(psth500k_b_single_fibers);
    
    for mc = 1:nrep(3)
       spont = sponts_concat(mc+nrep(1)+nrep(2));
       tabs = tabss_concat(mc+nrep(1)+nrep(2));
       trel = trels_concat(mc+nrep(1)+nrep(2));
       [psth500k_c_temp,~,~,~, ~,~] = model_Synapse_BEZ2018(vihc_temp,psth_freq(i),1,ts,1,0,spont,tabs,trel);
       psth500k_c_single_fibers(mc,:) = psth500k_c_temp;
    end

    psth500k_c = sum(psth500k_c_single_fibers);
    psth500k = psth500k_a + psth500k_b + psth500k_c;
    clear psth500k_a psth500k_b psth500k_c

    pr = sum(reshape(psth500k,binw,len),1)/sum(nrep)/binwidth; % psth in units of spikes/s/fiber
    pr = filtfilt(B,A,pr);
    psth_out = single(pr);
end
