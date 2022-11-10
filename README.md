---
title: "Non-invasive eye movement monitoring with eyes shut    "

---
[Slides](http://zijingzhang1997.github.io/files/eye_intro.pdf)

<img src='/images/ROG/ROG.gif'>  <br/>

This work presents radiooculogram (ROG), a novel sensor for non-invasive eye movement (EM) monitoring with eyes closed. We have experimentally demonstrated accurate measurements of EM frequency and directions for 5 participants and benchmarked ROG with electrooculogram (EOG). Compared with biopotential-based sensors, ROG has higher user comfort due to touchless operation and can capture direct muscle activity even in deep tissues. This work on voluntary EM sensing can serve as the baseline implementation for eventual sleep rapid EM monitoring.  

<img src='/images/ROG/pic1.png'>  <br/>
Experimental setup of ROG

<img src='/images/ROG/waveform.png'>  <br/>

<img src='/images/ROG/pic2.png'>  <br/>
EM frequency estimation

<img src='/images/ROG/pic4.png'>  <br/>
<img src='/images/ROG/pic3.png'>  <br/>
EM direction estimation

**Advantages**: 

•	Improved user comfort. EOG and EEG measurements demand numerous electrodes around the eye region with stable electrical contact, which are inconvenient, uncomfortable, and prone to face motion interference. ROG can operate without direct skin contact. For example, ROG can be integrated into the eyeglasses frame or implemented as an insert or at an outer layer of a sleep mask or a masquerade. As no skin contact electrode is involved and ROG sensors just need to be at a constant place relative to eyes, both comfort and design freedom are enhanced.

•	Unmediated sensing of directional EM. While the biopotential-based sensors such as EOG, EEG, and electromyography (EMG) measure neural stimulus for muscle activity, ROG directly measures the muscle motion by coupling RF energy to deep internal muscles, i.e., muscle condition is the direct output, instead of an estimate derived from the measured stimulus. EOG and ROG can also be used together to retrieve the closed loop of stimulation and actuation. 

•	Baseline for sleep REM detection. While camera-based methods are difficult to use for sleep REM, ROG has the flexibility to operate when eyes are open or closed without privacy concern. This work on conscious participants performing voluntary EM can formulate a validation baseline for future sleep REM monitoring.    

  
