
# Beam shaping using random phase-only mask for MATLAB®

* **Authors**: [M. Yommi](https://github.com/maxiyommi), [A. Bianchetti](https://github.com/abianchetti), [P. Etchepareborda](https://github.com/PablitoE) and [A. Federico](https://github.com/alefederico).
* **Email**: grupolatof@gmail.com
* **Homepage**: http://www.inti.gov.ar

## Reference
M. Yommi, A. Bianchetti, P. Etchepareborda, and A. Federico. Determining high-accuracy random phase-only masks for 
local complex modulation of light fields. 2019

> Please cite this paper when the provided code is used. See licensing terms for details.

## Contents
| Script  | Description |
|---|---|
| [main.m](/source/main.m) | Main script. |
| [define_T.m](/source/define_T.m.m) | Initial conditions settings. |
| [FFT2.m](/source/FFT.m) | Fast Fourier Transform balanced. |
| [iFFT2](/source/iFFT2.m) | Inverse Fast Fourier Transform balanced. |
| [PR.m](/source/PR.m) | Phase retrieval algorithm.  |
| [speckle_gen.m](/source/speckle_gen.m) |  Speckle field generator. |
| [ssim_index.m](/source/ssim_index.m)  |  Structural SIMilarity (SSIM) index. |
| [format_subplot.m](/source/format_subplot.m) | Subplot function.  |

For a detailed description of arguments and outputs consult the docstring in the files.
 
If you have suggestions, bugs or feature requests or want to contribute code, please email us.

## Disclaimer
All functions in this toolbox were implemented with care and tested on the examples presented in *Determining high-accuracy random phase-only masks for local complex modulation of light fields.* were possible. 
Nevertheless, they may contain errors or bugs, which may affect the outcome of your analysis. 
We do not take responsibility for any harm coming from using this toolbox, neither if it is caused by errors in the software nor if it is caused by its improper application. Please email us any bugs you find.

> Distributed under Open Source Script.




