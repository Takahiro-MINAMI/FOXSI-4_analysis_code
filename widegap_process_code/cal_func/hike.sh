#!/bin/zsh

## chmod 775 ./hike.sh

## record the main directory to work in
SCRIPT_DIR=$( cd -- "$( dirname -- "${(%):-%N}" )" &> /dev/null && pwd )

## split up the original root file path and extension
directory=$(dirname "$1")
file=$(basename "$1")
ext=$(echo $file | cut -d"." -f2)
filename=$(echo $file | cut -d"." -f1)

## define common files to use
basefilename=$directory"/"$filename
logfile=$basefilename"_hikeTerminalOutput.log"

function normal_line_gap() {
    ## a useful gap between the analysis chunks
    echo "\n\n" >> $logfile 2>&1
}

function normal_line() {
    ## save normal lines to the output file and run them
    echo "COMMAND: $@" >> $logfile 2>&1
    $@ >> $logfile 2>&1
}

function get_cdte_id() {
    ## CdTe ID mapping is
    ## **det_number {1, 2, 3, 4} corresponds to {"no2022_01", "no2021_07", "no2021_05", "no2021_06"}**
    if [[ "$1" == "1" ]]
    then
        cdteid="no2022_01"
    elif [[ "$1" == "2" ]]
    then
        cdteid="no2021_07"
    elif [[ "$1" == "3" ]]
    then
        cdteid="no2021_05"
    elif [[ "$1" == "4" ]]
    then
        cdteid="no2021_06"
    else
        echo "Unrecognised <cdte#>." >&2
        echo "Value should be in [1,4]." >&2
        exit 1
    fi
}

## define some default files to call
get_cdte_id $2
caldatatreefile=$SCRIPT_DIR/cal_func/"$cdteid"/enCalDataTree.root
thresholdfile=$SCRIPT_DIR/cal_func/"$cdteid"/threth_keV.txt
respfile=$SCRIPT_DIR"/gaploss_resp/resp_"$cdteid"_low.root"
doirespfile=$SCRIPT_DIR/doiloss_resp/doi_resp_"$cdteid".root

## some code from https://medium.com/@wujido20/handling-flags-in-bash-scripts-4b06b4d0ed04
# Function to display script usage
function usage() {
    echo ""
    echo "Usage: $0 [OPTIONS]"
    echo ""
    echo "path/to/hike.sh <path_and_name_of_raw_file> <cdte#> [OPTIONS]"
    echo ""
    echo "Options:"
    echo " -h,  --help                 Display this help message"
    echo " -c*, --caldatatreefile*     Define calibration file"
    echo " -t*, --thresholdfile*       Define threshold file"
    echo " -r*, --respfile*            Define response file"
    echo " -d*, --doirespfile*         Define DOI file"
    echo ""
    echo "Generic Eamples:"
    echo " path/to/hike.sh <path_and_name_of_raw_file> <cdte#> -h"
    echo " path/to/hike.sh <path_and_name_of_raw_file> <cdte#> -c <cal_file>"
    echo " path/to/hike.sh <path_and_name_of_raw_file> <cdte#> -c=<cal_file>"
    echo " path/to/hike.sh <path_and_name_of_raw_file> <cdte#> --caldatatreefile <cal_file>"
    echo " path/to/hike.sh <path_and_name_of_raw_file> <cdte#> --caldatatreefile=<cal_file>"
    echo ""
    echo "Other Eamples:"
    echo " ## For CdTe1:"
    echo " path/to/hike.sh path/to/raw/root/file.root 1 -c path/to/cal/file.root"
}

function has_argument() {
    ## check if the flag has a value
    [[ ("$1" == *=* && -n ${1#*=}) || ( ! -z "$2" && "$2" != -*)  ]];
}

function extract_argument() {
    ## get the value from a flag
    echo "${2:-${1#*=}}"
}

# Function to handle options and arguments
function handle_options() {
    ## handles the optional inputs to the script
    while [ $# -gt 0 ]; do
        case $1 in
            -h | --help)
                usage
                exit 0
                ;;
            -c* | --caldatatreefile*)
                if ! has_argument $@; then
                    echo "Calibration file not specified." >&2
                    usage
                    exit 1
                fi
                caldatatreefile=$(extract_argument $@)
                ;;
            -t* | --thresholdfile*)
                if ! has_argument $@; then
                    echo "Threshold file not specified." >&2
                    usage
                    exit 1
                fi
                thresholdfile=$(extract_argument $@)
                ;;
            -r* | --respfile*)
                if ! has_argument $@; then
                    echo "Response file not specified." >&2
                    usage
                    exit 1
                fi
                respfile=$(extract_argument $@)
                ;;
            -d* | --doirespfile*)
                if ! has_argument $@; then
                    echo "DOI file not specified." >&2
                    usage
                    exit 1
                fi
                doirespfile=$(extract_argument $@)
                ;;
        esac
        shift
    done
}

## go through the inputs and set things up
handle_options "$@"

echo "Let's go on a hike!" | tee $logfile 2>&1
## record the date and time
echo "Date & Time: "$(date +%d.%m.%y-%H:%M:%S) | tee -a $logfile 2>&1
## need to change directory because of relative paths being called
normal_line cd $SCRIPT_DIR
normal_line_gap

## make sure to log the files being called into the pipeline
echo "Files being called are:" >> $logfile 2>&1
echo "caldatatreefile=$caldatatreefile" >> $logfile 2>&1
echo "thresholdfile=$thresholdfile" >> $logfile 2>&1
echo "respfile=$respfile" >> $logfile 2>&1
echo "doirespfile=$doirespfile" >> $logfile 2>&1
normal_line_gap

## set-up is done, now actually run the PATH scripts
echo "PATH 1" | tee -a $logfile 2>&1
echo $SCRIPT_DIR/path1_base_forwidegap $1 $2
normal_line $SCRIPT_DIR/path1_base_forwidegap $1 $2
normal_line_gap
echo "PATH 2" | tee -a $logfile 2>&1
normal_line $SCRIPT_DIR/path2_cal_and_merge_eachside $basefilename"_base.root" $caldatatreefile $thresholdfile
normal_line_gap
echo "PATH 3" | tee -a $logfile 2>&1
normal_line $SCRIPT_DIR/path3_mod_gaploss $basefilename"_base_cal.root" $respfile
normal_line_gap
echo "PATH 4" | tee -a $logfile 2>&1
normal_line $SCRIPT_DIR/path4_mod_doiloss $basefilename"_base_cal_gapmod.root" $doirespfile
normal_line_gap

## want to move the files into a directory structure
echo "Moving outputs into directory structure located in:" | tee -a $logfile 2>&1
echo "    $directory" | tee -a $logfile 2>&1
echo "" >> $logfile 2>&1
MOVING_STUFF_LINE="Moving stuff around, beep beep."
## path 1
echo $MOVING_STUFF_LINE": PATH 1" >> $logfile 2>&1
normal_line mkdir $directory"/path1"
normal_line mv $basefilename"_base.root" $directory"/path1/"$filename"_base.root"
normal_line_gap
## path 2
echo $MOVING_STUFF_LINE": PATH 2" >> $logfile 2>&1
normal_line mkdir $directory"/path2"
normal_line mv $basefilename"_base_cal.root" $directory"/path2/"$filename"_base_cal.root"
normal_line mkdir $directory"/path2/figs"
normal_line mv $basefilename"_base_cal.pdf" $directory"/path2/figs/"$filename"_base_cal.pdf"
normal_line_gap
## path 3
echo $MOVING_STUFF_LINE": PATH 3" >> $logfile 2>&1
normal_line mkdir $directory"/path3"
normal_line mv $basefilename"_base_cal_gapmod.root" $directory"/path3/"$filename"_base_cal_gapmod.root"
normal_line mkdir $directory"/path3/figs"
normal_line mv $basefilename"_base_cal_gapmod.pdf" $directory"/path3/figs/"$filename"_base_cal_gapmod.pdf"
normal_line_gap
## path 4
echo $MOVING_STUFF_LINE": PATH 4" >> $logfile 2>&1
normal_line mkdir $directory"/path4"
normal_line mv $basefilename"_base_cal_gapmod_doimod.root" $directory"/path4/"$filename"_base_cal_gapmod_doimod.root"
normal_line mkdir $directory"/path4/figs"
normal_line mv $basefilename"_base_cal_gapmod_doimod.pdf" $directory"/path4/figs/"$filename"_base_cal_gapmod_doimod.pdf"
normal_line_gap

echo $MOVING_STUFF_LINE": PATH 5" >> $logfile 2>&1
normal_line mkdir $directory"/path5"
normal_line mv $basefilename"_base_cal_gapmod_doimod_poscal.root" $directory"/path5/"$filename"_base_cal_gapmod_doimod"_poscal.root"
normal_line mkdir $directory"/path5/figs"
normal_line mv $basefilename"_base_cal_gapmod_doimod_poscal.pdf" $directory"/path5/figs/"$filename"_base_cal_gapmod_doimod_poscal.pdf"
normal_line_gap

echo "Hike finished. Now get sciencing if all is OK in the world!" | tee -a $logfile 2>&1
