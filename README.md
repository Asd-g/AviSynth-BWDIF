# Description

Motion adaptive deinterlacing based on yadif with the use of w3fdif and cubic interpolation algorithms.

This is [a port of the VapourSynth plugin Bwdif](https://github.com/HomeOfVapourSynthEvolution/VapourSynth-Bwdif).

# Usage

```
BWDIF (clip, int "field", clip "edeint", int "opt")
```

## Parameters:

- clip\
    A clip to process. All planar formats are supported.
    
- field\
    Controls the mode of operation (double vs same rate) and which field is kept.\
    -2: Double rate (alternates each frame), AviSynth internal order.\
    -1: Same rate, AviSynth internal order.\
    0: Same rate, keep bottom field.\
    1: Same rate, keep top field.\
    2: Double rate (alternates each frame), starts with bottom.\
    3: Double rate (alternates each frame), starts with top.\
    Default: -1.
    
- edeint\
    Clip from which to take spatial predictions. This clip must be the same width, height, and colorspace as the input clip.\
    If using same rate output, this clip should have the same number of frames as the input. If using double rate output, this clip should have twice as many frames as the input.
    
- opt\
    Sets which cpu optimizations to use.\
    -1: Auto-detect.\
    0: Use C++ code.\
    1: Use SSE2 code.\
    2: Use AVX2 code.\
    3: Use AVX512 code.\
    Default: -1.