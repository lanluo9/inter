for imouse = 1:3
    resp_ad(imouse) = max(trace_mouse_avg(imouse, 1:20)) - min(trace_mouse_avg(imouse, 1:20));
    resp_ad_targ(imouse) = max(trace_mouse_avg(imouse, 21:60)) - min(trace_mouse_avg(imouse, 21:60));
    adp_mouse_trace(imouse) = resp_ad_targ(imouse)/resp_ad(imouse) - 1;
end
adp_mouse_trace