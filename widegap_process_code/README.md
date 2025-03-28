# Process code for FOXSI-4 wide-gap CdTe-DSD

## Reminded Tasks:  
- Implementation for sub-strip position determination
- Fine position correction using Pointing information 
- Handling multiplicity events
- Low-energy threshold for different detectors
- Gain-shifts (may caused by canister temp diff.?)  


## Updated Memo:  
##### 2024.08.18.  
- Cal/Gaploss/DoI functions for each flight detector are updated.  
- Bug related to the primary threshold in path2 is modified.
- Event quality flag is added after path3 for handling multiplicity events.  
  
##### 2024.06.30.  
- DoI cal functions for each flight detector are added.  
- Calibration functions for each flight detector are updated (adjust for Am241 60 keV peak).  
 
##### 2024.06.24.  
- Programs for making cal function (enCalDataTree.root) are added.  
  
##### 2024.06.16.  
- Position calibration in solar coordinates (helio_*_arcsec) is implemented in path5 process based on WSMR test  
- Position cal files (detector coordinate -> arcsec) are added in poscal_resp folder  
- Threshold for charge-share is changed to 2.0 keV considering the noise level (need to be discussed)  
  
##### 2024.06.01.  
- PPS-based time synchronization is implemented in path1 process
- Time cal files (ti->PPS, PPS->Realtime) are added in cal_func folder
- Add Makefile for easy compile  
  
##### 2024.05.27.  
- Initial version  
  

## Overview:
These are source codes for FOXSI-4 wide-gap CdTe-DSDs.    
Each step of the code is written in C++. Therefore, please compile it when used.  
ex)   
```
clang++ -O3 -mtune=native -march=native  `root-config --libs --cflags` -Wl,-rpath,/Users/shunsaku/soft/root/root-6_30_06_install/lib/ -o path1_base_forwidegap path1_base_forwidegap.cpp
```
or just type 
```
make
```
The path to ROOT should be referred as ${ROOTSYS}.  
det_number {1, 2, 3, 4} corresponds to {"no2022_01", "no2021_07", "no2021_05", "no2021_06"}.  

#### STEP0: rawdata2root
-> test.root  

#### STEP1: path1_base_forwidegap.cpp
##### Remap ch, CMN subtraction, Time Sync.  
Input: test.root  
Output: test_base.root  
USAGE:   
./base_forwidegap 'input filename' 'det_number' flag_all  
ex)   
```
./base_forwidegap ./test.root #Non Flight data  
./base_forwidegap ./test.root 1 #For CdTe #1 Flight data
./base_forwidegap ./test.root 1 -all  
```
#### STEP2: path2_cal_and_merge_eachside.cpp
##### Calibration for each ch/side (adc->Energy).  
Input: test_base.root, cal-func. rootfile, threth. txtfile  
Output: test_base_cal.root  
USAGE:  
./path2_cal_and_merge_eachside 'input filename' 'cal func. filename' 'threshold filename' -all  
ex)   
```
./cal_and_merge_eachside test_base.root cal_func/no2021_05/enCalDataTree.root cal_func/no2021_05/threth_keV.txt  
```

#### STEP3: path3_mod_gaploss.cpp
##### Correction for wide-gap loss for each side  
Input: test_base_cal.root, gaploss-func. rootfile  
Output: test_base_cal_gapmod.root  
USAGE:  
./path3_mod_gaploss 'input filename' 'gaploss func. filename' -all  
ex)   
```
./path3_mod_gaploss test_base_cal.root gaploss_resp/resp_no2021_05.root  
```

#### STEP4: path4_mod_doiloss.cpp  
##### Correction for DoI loss using both sides of energies.  
##### CAUTION: If there is a gain-shift, this method will not work correctly.
Input: test_base_cal_gapmod.root, doiloss-func. rootfile  
Output: test_base_cal_gapmod_doiloss.root  
USAGE:   
./path4_mod_doiloss 'input filename' 'doiloss func. filename'  
ex)   
```
./path4_mod_doiloss test_base_cal_gapmod.root doiloss_resp/doi_resp_no2021_02.root  
```

#### STEP5: path5_poscal.cpp  
##### Position calibration in solar coordinates (helio_*_arcsec) 
Input: test_base_cal_gapmod_doiloss.root, poscal-func. rootfile  
Output: test_base_cal_gapmod_doiloss_poscal.root  
USAGE:   
./path5_poscal 'input filename' 'poscal func. filename' 'det_number'
ex)   
```
./path5_poscal test_base_cal_gapmod_doiloss.root poscal_resp/position_20240514.root 1  
```

#### Running `./hike.sh`

This script is to help run all the above paths for a given raw file rather than runnig all the path steps individually. The `hike.sh` script should be runnable from any directory.

The script will have to be executable which means the user will like have to do _one_ of the following (or somthing similar):

- `chmod 775 hike.sh`
- `chmod +x hike.sh`

The usage is as follows (can can be found by running`path/to/hike.sh <path_and_name_of_raw_file> <cdte#> -h` where `<cdte#>` is an integer in `[1,4]`):

```c
path/to/hike.sh <path_and_name_of_raw_file> <cdte#> [OPTIONS]

Options:
 -h,  --help                 Display this help message"
 -c*, --caldatatreefile*     Define calibration file"
 -t*, --thresholdfile*       Define threshold file"
 -r*, --respfile*            Define response file"
 -d*, --doirespfile*         Define DOI file"

Generic Eamples:
 path/to/hike.sh <path_and_name_of_raw_file> <cdte#> -h
 path/to/hike.sh <path_and_name_of_raw_file> <cdte#> -c <cal_file>
 path/to/hike.sh <path_and_name_of_raw_file> <cdte#> -c=<cal_file>
 path/to/hike.sh <path_and_name_of_raw_file> <cdte#> --caldatatreefile <cal_file>
 path/to/hike.sh <path_and_name_of_raw_file> <cdte#> --caldatatreefile=<cal_file>

Other Eamples:
 ## For CdTe1:
 path/to/hike.sh path/to/raw/root/file.root 1 -c path/to/cal/file.root
```
