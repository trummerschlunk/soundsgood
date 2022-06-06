// double precision -double needed!
//
import("stdfaust.lib");

// init values

Nch = 2; //number of channels

init_noisegate_threshold = -70; // not used in voc version

init_leveler_target = -16;
init_leveler_maxboost = 40;
init_leveler_maxcut = 40;
init_leveler_gatethreshold = -45;
init_leveler_speed = .05;

init_mscomp10_postgain = 1;

init_kneecomp_thresh = 2;
init_kneecomp_postgain = 0;

init_limiter_lad_ceil = -2;
init_limiter_postgain = 0;

init_brickwall_ceiling = -1;


target = hslider("../../[1]TARGET[unit:dB]", init_leveler_target,-50,0,1);


// main
process =

    // ba.bypass2(checkbox("bypass all"),

    hgroup("row1",

    si.bus(2) :

    //hgroup("MASTER_ME", hgroup("[0]INPUT",peak_meter(Nch))) :
    //hgroup("MASTER_ME", hgroup("[0]INPUT",lufs_any(Nch))) :
    //hgroup("MASTER_ME", hgroup("[1]STEREO CORRECT",correlate_meter)) :
    //hgroup("[1]STEREO CORRECT",correlate_correct_bp) :



    dc_filter(2) :

    hgroup("[2]NOISEGATE",noisegate(2)):

    hgroup("[3]LEVELER",leveler(target))


    // hgroup("MASTER_ME", vgroup("[3]MULTIBAND COMPRESSOR", mbcomp(Nch))) :

    ) :  // end vgroup row1

    hgroup("row2",
    hgroup("[3]MSCOMP10", mscomp10(target))  :

    // hgroup("row3",
    // hgroup("[3]5-BAND COMPRESSOR", comp5st) :
    // hgroup("[3]MULTIBAND MID-SIDE COMPRESSOR", mbmscomp(Nch)) :
    hgroup("[4]KNEECOMP",kneecomp(target)) :

    hgroup("[7]LIMITER", limiter) :
    hgroup("[8]BRICKWALL",brickwall) :
    hgroup("[9]OUTPUT",lufs_any(Nch)) :
    //hgroup("MASTER_ME", hgroup("[9]OUTPUT",peak_meter(Nch))) :


    si.bus(2)
    ) // end vgroup row2
    //) // hgroup MASTER_Me end
    // ) // bypass end
;


// DC FILTER
dc_filter(N) = par(i,N,fi.dcblocker);

// NOISE GATE
noisegate(N) = gate_any(N,noisegate_thresh,noisegate_attack,noisegate_hold,noisegate_release) with {
    noisegate_thresh = vslider("[0]threshold",init_noisegate_threshold, -95, 0, 1);
    noisegate_attack = 0.01;
    noisegate_hold = 1;
    noisegate_release = 2;

    gate_any(N,thresh,att,hold,rel) = B <: B, (B :> ggm : vbargraph("[2]gate level",0,1) <: B) : ro.interleave(N,2) : par(i,N,*)
    with {
        B = si.bus(N);
        ggm = gate_gain_mono(thresh,att,hold,rel);
    };

    gate_gain_mono(thresh,att,hold,rel,x) = x : extendedrawgate : an.amp_follower_ar(att,rel) with {
        extendedrawgate(x) = max(float(rawgatesig(x)),holdsig(x));
        rawgatesig(x) = inlevel(x) > ba.db2linear(thresh);
        minrate = min(att,rel);
        inlevel = an.amp_follower_ar(minrate,minrate);
        holdcounter(x) = (max(holdreset(x) * holdsamps,_) ~-(1));
        holdsig(x) = holdcounter(x) > 0;
        holdreset(x) = rawgatesig(x) < rawgatesig(x)'; // reset hold when raw gate falls
        holdsamps = int(hold*ma.SR);
    };
};



// LEVELER
leveler(target) = B <:    B    ,   (B <: B,B : LkN, + : calc : _ <: B) : ro.interleave(N,2) : par(i,N,*)
    with {
    N = 2;
    B = si.bus(N);

    calc(lufs,sc) = (lufs : vbargraph("[1][unit:dB]LUFS",-70,0) : (target - _) : lp1p(leveler_speed_gated(sc)) : limit(limit_neg,limit_pos) : vbargraph("[2]gain",-50,50) : ba.db2linear) , sc : _,!;


    // target = vslider("[3]target[unit:dB]", init_leveler_target,-50,0,1);

    limit_pos = vslider("[5]max +", init_leveler_maxboost, 0, 60, 1);
    limit_neg = vslider("[6]max -", init_leveler_maxcut, 0, 60, 1) : ma.neg;
    limit(lo,hi) = min(hi) : max(lo);

    leveler_speed = vslider("[4]speed", init_leveler_speed, .005, 0.15, .005);
    leveler_speed_gated(sc) = (gate_gain_mono(leveler_gate_thresh,0.1,0,0.1,abs(sc)) <: attach(_, (1-_) : vbargraph("[7]leveler gate",0,1))) : _ * leveler_speed;
    leveler_gate_thresh = vslider("[8]lev gate thresh[unit:dB]", init_leveler_gatethreshold,-90,0,1);

    // from library:
    gate_gain_mono(thresh,att,hold,rel,x) = x : extendedrawgate : an.amp_follower_ar(att,rel) with {
        extendedrawgate(x) = max(float(rawgatesig(x)),holdsig(x));
        rawgatesig(x) = inlevel(x) > ba.db2linear(thresh);
        minrate = min(att,rel);
        inlevel = an.amp_follower_ar(minrate,minrate);
        holdcounter(x) = (max(holdreset(x) * holdsamps,_) ~-(1));
        holdsig(x) = holdcounter(x) > 0;
        holdreset(x) = rawgatesig(x) < rawgatesig(x)'; // reset hold when raw gate falls
        holdsamps = int(hold*ma.SR);
    };

};



// 10BAND MID-SIDE COMPRESSOR
mscomp10(target) = _,_ : ms_enc : par(i,2,fibank_mono) : ro.interleave(N,2) : par(i,N,compst(rdtable(thresh_offset,i))) : par(i,N,ms_dec) :> _,_ : post_gain with{

    // threshold offset (high freq to low freq)
    thresh_offset = waveform{-25,-17,-14,-13,-11,-10,-8,0,0,-6};

    M = 1;
    ftop = 10000;
    N = 10 * M;
    fibank_mono = fi.mth_octave_filterbank_default(M,ftop,N);

    // stereo to m/s encoder
    ms_enc = _*0.5,_*0.5 <: +, -;
    // m/s to stereo decoder
    ms_dec = _,_ <: +, -;

    // stereo compressor
    compst(thr_os) = co.FBcompressor_N_chan(strength,thresh+thr_os,att,rel,knee,prePost,link,meter,2) with {
        strength = 0.1;
        thresh = target + vslider("[unit:dB]tar-thr",-2,-10,10,1);
        att = 0.015;
        rel = 0.6;
        knee = 12;
        prePost = 1;
        link = 0.5;
        meter = _; //_ <: (_, (ba.linear2db : ma.neg : vbargraph("[unit:dB]",0,3))) : attach;
    };

    //post_gain
    post_gain = par(i,2,_ * g) with {
        g =  vslider("post gain[unit:dB]", init_mscomp10_postgain,-10,+10,0.5) : ba.db2linear;
    };
};




// KNEE COMPRESSOR

kneecomp(target) = ms_enc : co.RMS_FBcompressor_peak_limiter_N_chan(strength,thresh,threshLim,att,rel,knee,link,meter,meterLim,2) : ms_dec : post_gain with {

  strength = 0.1; //vslider("strength", 0.1, 0, 1, 0.1);
  thresh = target + vslider("[unit:dB]tar-thr",init_kneecomp_thresh,-12,6,1);
  threshLim = +3; //vslider("threshLim",3,-12,3,1);
  att = 0.4; //vslider("att",0.4,0.001,1,0.001);
  rel = 0.8; //vslider("rel",0.8,0.01,1,0.001);
  knee = 12; //vslider("knee",12,0,12,1);
  link = 0.5; //vslider("link", 0.5, 0, 1, 0.1);
  meter = _<: _,(ba.linear2db : ma.neg : vbargraph("[unit:dB]",0,3)) : attach;
  meterLim = _<: _,(ba.linear2db : ma.neg : vbargraph("[unit:dB]",0,3)) : attach;

  // stereo to m/s encoder
  ms_enc = _*0.5,_*0.5 <: +, -;
  // m/s to stereo decoder
  ms_dec = _,_ <: +, -;

  //post_gain
    post_gain = par(i,2,_ * g) with {
        g =  vslider("post gain[unit:dB]", init_kneecomp_postgain,-10,+10,0.5) : ba.db2linear;
    };

};







// LIMITER
limiter = limiter_lad_N(2,limiter_lad_lookahead, init_limiter_lad_ceil : ba.db2linear, limiter_lad_attack, limiter_lad_hold, limiter_lad_release) : post_gain with{

    N=2;

    limiter_lad_lookahead = 0;
    limiter_lad_attack = 0.001;
    limiter_lad_hold = 0.005;
    limiter_lad_release = 0.02;

    // lookahead limiter (N-channel)
    limiter_lad_N(N, LD, ceiling, attack, hold, release) =
        si.bus(N) <: par(i, N, @ (LD * ma.SR)),
            (scaling <: si.bus(N)) : ro.interleave(N, 2) : par(i, N, *)
        with {
            scaling = ceiling / max(amp_profile, ma.EPSILON) : min(1) : meter_limiter_lad_N;
            amp_profile = par(i, N, abs) : maxN(N) : ba.peakholder(hold * ma.SR) : att_smooth(attack) : rel_smooth(release);
            att_smooth(time, in) = si.smooth(ba.tau2pole(time), in);
            rel_smooth(time, in) = an.peak_envelope(time, in);
            maxN(1) = _;
            maxN(2) = max;
            maxN(N) = max(maxN(N - 1));
        };

    // post_gain
    post_gain = par(i,Nch,_ * g) with {
        g =  vslider("[9]post gain[unit:dB]", init_limiter_postgain,-10,+10,0.5) : ba.db2linear;
    };

    // metering
    //meter_limiter_lad_N = _ <: attach(ba.linear2db : vbargraph("[8][unit:dB]GR",-12,0));
    meter_limiter_lad_N = _ <: attach(ba.linear2db : abs : vbargraph("[8][unit:dB]GR",0,12));
};


// BRICKWALL
brickwall = limiter_lad_N(N, limiter_lad_lookahead, limiter_lad_ceil, limiter_lad_attack, limiter_lad_hold, limiter_lad_release)
    with{

    N=2;

    twopi = 2 * ma.PI;

    limiter_lad_lookahead = 0.01;
    limiter_lad_ceil = init_brickwall_ceiling : ba.db2linear;
    limiter_lad_attack = .01 / twopi;
    limiter_lad_hold = .1;
    limiter_lad_release = 1 / twopi;

    // lookahead limiter (N-channel)
    limiter_lad_N(N, LD, ceiling, attack, hold, release) =
        si.bus(N) <: par(i, N, @ (LD * ma.SR)),
            (scaling <: si.bus(N)) : ro.interleave(N, 2) : par(i, N, *)
        with {
            scaling = ceiling / max(amp_profile, ma.EPSILON) : min(1) : meter_limiter_lad_N;
            amp_profile = par(i, N, abs) : maxN(N) : ba.peakholder(hold * ma.SR) : att_smooth(attack) : rel_smooth(release);
            att_smooth(time, in) = si.smooth(ba.tau2pole(time), in);
            rel_smooth(time, in) = an.peak_envelope(time, in);
            maxN(1) = _;
            maxN(2) = max;
            maxN(N) = max(maxN(N - 1));
        };

    // metering
    meter_limiter_lad_N = _ <: attach(ba.linear2db : abs : vbargraph("[8][unit:dB]GR",0,12));
};

// METERING
peak_meter(N) = par(i, N, (_ <: attach(_, envelop : vbargraph("[unit:dB]CH%i", -70, 0)))) with{

    vmeter(x) = attach(x, envelop(x) : vbargraph("[unit:dB]", -70, 0));
    envelop = abs : max(ba.db2linear(-70)) : ba.linear2db : min(10) : max ~ -(40.0/ma.SR);
};

// +++++++++++++++++++++++++ LUFS METER +++++++++++++++++++++++++
Tg = 3; // 3 second window for 'short-term' measurement
zi = an.ms_envelope_rect(Tg); // mean square: average power = energy/Tg = integral of squared signal / Tg
kfilter = fi.highpass(1, 60) : fi.high_shelf(4, 1800);


// 2-channel
lk2 = par(i,2,kfilter : zi) :> 10 * log10(max(ma.EPSILON)) : -(0.691);


//envelope via lp by Dario Sanphilippo
lp1p(cf, x) = fi.pole(b, x * (1 - b)) with {
    b = exp(-2 * ma.PI * cf / ma.SR);
};
// zi_lp(x) = lp1p(1 / Tg, x * x);

// one channel
Lk = kfilter : zi : 10 * log10(max(ma.EPSILON)) : -(0.691);

// N-channel
LkN = par(i,Nch, kfilter : zi) :> 10 * log10(max(ma.EPSILON)) : -(0.691);

// N-channel by Yann Orlarey
lufs_any(N) = B <: B, (B : par(i,N,kfilter:zi) :> 10 * log10(max(ma.EPSILON)) : -(0.691) : vbargraph("[unit:dB]LUFS",-70,0)) : si.bus(N-1), attach(_,_)
    with {
        B = si.bus(N);

    };


// correlation meter
correlate_meter(x,y) = x,y <: x , attach(y, (corr(t) : vbargraph("correlation",-1,1))) : _,_ with {
    t = .2; // averaging period in seconds

    avg(t, x) = fi.pole(p, (1 - p) * x) // 1-pole lowpass as average
        with {
            p = exp((((-2.0 * ma.PI) / t) / ma.SR));
        };
    var(t, x) = avg(t, (x - avg(t, x)) ^ 2); // variance
    sd(t, x) = sqrt(var(t, x)); // standard deviation
    cov(t, x1, x2) = avg(t, (x1 - avg(t, x1)) * (x2 - avg(t, x2))); // covariance
    corr(t, x1, x2) = cov(t, x1, x2) / max(ma.EPSILON, (sd(t, x1) * sd(t, x2))); // correlation
};

// stereo correction based on correlation
correlate_correct(l,r) = out_pos1, out_neg1, out_0, out_pos, out_neg :> _,_ with {

    t = .2; // averaging period in seconds

    avg(t, x) = fi.pole(p, (1 - p) * x) // 1-pole lowpass as average
        with {
            p = exp((((-2.0 * ma.PI) / t) / ma.SR));
        };
    var(t, x) = avg(t, (x - avg(t, x)) ^ 2); // variance
    sd(t, x) = sqrt(var(t, x)); // standard deviation
    cov(t, x1, x2) = avg(t, (x1 - avg(t, x1)) * (x2 - avg(t, x2))); // covariance
    corr(t, x1, x2) = cov(t, x1, x2) / max(ma.EPSILON, (sd(t, x1) * sd(t, x2))); // correlation


    th =.0001;
    corr_pos1 = avg(t, (corr(t,l,r) > (1-th))) : smoothing /*: vbargraph("[5]1",0,1)*/;
    corr_neg1 = avg(t, corr(t,l,r) < (-1+th)) : smoothing /*: vbargraph("[9]-1",0,1)*/;
    corr_0 = avg(t, ((corr(t,l,r) < th) & (corr(t,l,r) > (0-th)))) : smoothing /*: vbargraph("[7]0",0,1)*/;
    corr_pos = avg(t, ((corr(t,l,r) > (0+th)) & (corr(t,l,r) < (1-th)))) : smoothing /*: vbargraph("[6]>0,<1",0,1)*/;
    corr_neg = avg(t, ((corr(t,l,r) > (-1+th)) & (corr(t,l,r) < (0-th)))) : smoothing /*: vbargraph("[8]>-1,<0",0,1)*/;

    smoothing = lp1p(2) ;

    out_pos1 = ((l * corr_pos1 + r * corr_pos1) /2) , ((l * corr_pos1 + r * corr_pos1) /2);
    out_neg1 = ((l * corr_neg1 + (-r) * corr_neg1) /2) , ((l * corr_neg1 + (-r) * corr_neg1) /2);
    out_0 = (l * corr_0 + r * corr_0) , (l * corr_0 + r * corr_0);
    out_pos = l * corr_pos , r * corr_pos;
    out_neg = l * corr_neg , r * corr_neg; // old: out_neg = l * corr_neg , (0-(r * corr_neg));
};

// stereo correction bypass checkbox
correlate_correct_bp = ba.bypass2(checkbox("bypass"), correlate_correct);
