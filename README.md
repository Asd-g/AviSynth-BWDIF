## Description

Motion adaptive deinterlacing based on yadif with the use of w3fdif and cubic interpolation algorithms.

This is [a port of the VapourSynth plugin Bwdif](https://github.com/HomeOfVapourSynthEvolution/VapourSynth-Bwdif).

### Requirements:

- AviSynth 2.60 / AviSynth+ 3.4 or later

- Microsoft VisualC++ Redistributable Package 2022 (can be downloaded from [here](https://github.com/abbodi1406/vcredist/releases)) (Windows only)

### Usage:

```
BWDIF (clip, int "field", clip "edeint", int "opt", float "thr", bool "debug", bool "pass")
```

### Parameters:

- clip\
    A clip to process. All planar formats are supported.

- field\
    Controls the mode of operation (double vs same rate) and which field is kept.\
    -2: Double rate (alternates each frame). If `_FieldBased` > `0` - `_FieldBased` frame property order. If ` _FieldBased` is missing or not supported - AviSynth internal order.\
    -1: Same rate. If `_FieldBased` > `0` - `_FieldBased` frame property order. If ` _FieldBased` is missing or not supported - AviSynth internal order.\
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

- thr\
    Threshold for interpolation.\
    If the difference between pixels of the prev/next frame is less than or equal to this, the resulted pixel wouldn't be interpolated.\
    Must be between 0.0..100.0.\
    100.0: No interpolation is performed.\
    Default: 0.0.

- debug\
    Whether to show which pixels will be interpolated.\
    Default: False.

- pass\
    Whether to return the source frame (repeated when double rate) when `_FieldBased` is `0`.
    Default: False.

### Building:

- Windows\
    Use solution files.

- Linux\
    ```
    Requirements:
        - Git
        - C++17 compiler
        - CMake >= 3.16
    ```
    ```
    git clone https://github.com/Asd-g/AviSynth-BWDIF && \
    cd AviSynth-BWDIF && \
    mkdir build && \
    cd build && \

    cmake ..
    make -j$(nproc)
    sudo make install
    ```
