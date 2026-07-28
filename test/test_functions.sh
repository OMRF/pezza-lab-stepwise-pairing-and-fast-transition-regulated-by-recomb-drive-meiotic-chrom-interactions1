#!/bin/bash

TEST_NUM=0
TEST_PASSED=0
SHORT_PATH=${SHORT_PATH:=FALSE}
DEFAULT_MIN_FSTRCMP_SCORE=${DEFAULT_MIN_FSTRCMP_SCORE:="0.97"}
DEBUG=${DEBUG:=FALSE}

failed () {
    RESULT=$1
    LEVELS=${2:-1} # defaults to 1 (function within which this is called)
    TEST_NUM=$((TEST_NUM+1))
    echo -e "not ok $TEST_NUM - (${FUNCNAME[$LEVELS]}) $RESULT"
}

passed () {
    RESULT=$1
    LEVELS=${2:-1} # defaults to 1 (function within which this is called)
    TEST_NUM=$((TEST_NUM+1))
    TEST_PASSED=$((TEST_PASSED+1))
    echo -e "ok $TEST_NUM - (${FUNCNAME[$LEVELS]}) $RESULT"
}

passed_file_is_as_expected () {
    result_file=$1
    LEVELS=${2:-2} # Defaults to function that called this one
    if [[ "$SHORT_PATH" == "TRUE" ]]; then
        passed "$( basename $result_file) is as expected." $LEVELS
    else
        passed "$result_file is as expected." $LEVELS
    fi
}

get_difference_between () {
    set +E # Turn off error trapping in functions (so we can discover a
           #   difference without choking on an error)
    DIFFERENCE=$(cmp $1 $2; exit 0)
    echo $DIFFERENCE
}

compare_ignoring_comments () {
    result_file=$1
    expected_file=$2 

    # First check if the file exists
    if [[ ! -s $result_file ]]; then
        failed "File empty or does not exist: $result_file"
    else
        no_comments_result=$(mktemp)
        grep -v '^#' $result_file > $no_comments_result
        no_comments_expected=$(mktemp)
        grep -v '^#' $expected_file > $no_comments_expected
        DIFFERENCE=$(get_difference_between $no_comments_result $no_comments_expected)
        if [[ -z "${DIFFERENCE}" ]]; then
            passed_file_is_as_expected "$result_file"
        else
            if [[ -s ${expected_file}.alt ]]; then
                compare_ignoring_comments $result_file ${expected_file}.alt # Yes recursion works in Bash!
            else
                failed "'${DIFFERENCE}'"
            fi
        fi
    fi
}


compare_temp_copies () {
    result_file=$1       # File to compare against expected
    expected_file=$2     # File to compare result against
    temp_result=$3
    temp_expected=$4

    # First check if the file exists
    if [[ ! -s $result_file ]]; then
        failed "File empty or does not exist: $result_file" 2
    else
        DIFFERENCE=$(get_difference_between $temp_result $temp_expected)
        if [[ -z "${DIFFERENCE}" ]]; then
            passed_file_is_as_expected "$result_file" 3 # Get level back to original calling function
        else
            if [[ -s ${expected_file}.alt ]]; then
                test_text_file $result_file ${expected_file}.alt # Yes recursion works in Bash!
            else
                failed "'$result_file' differed more than expected from '$expected_file': Approximated copies $DIFFERENCE" 2
            fi
        fi
    fi
}

test_text_file () {
    result_file=$1       # File to compare against expected
    expected_file=$2     # File to compare result against

    # First check if the file exists
    if [[ ! -s $result_file ]]; then
        failed "File empty or does not exist: $result_file"
    else
        DIFFERENCE=$(get_difference_between $result_file $expected_file)
        if [[ -z "${DIFFERENCE}" ]]; then
            passed_file_is_as_expected "$result_file"
        else
            if [[ -s ${expected_file}.alt ]]; then
                test_text_file $result_file ${expected_file}.alt # Yes recursion works in Bash!
            else
                failed "'${DIFFERENCE}'"
            fi
        fi
    fi
}

test_approx_sorted_text () {

    result_file=$1       # File to compare against expected
    expected_file=$2     # File to compare result against
    min_score=${3:-$DEFAULT_MIN_FSTRCMP_SCORE}
    
    # First check if the file exists
    if [[ ! -s $result_file ]]; then
        failed "File empty or does not exist: $result_file"
    else
        # Make temp files unreadable to outsiders
        umask 077
        
        # Create temp files
        result_sorted=$(mktemp)
        expected_sorted=$(mktemp)
        
        # Capture sorted versions of input files
        sort $result_file > $result_sorted
        sort $expected_file > $expected_sorted
        
        # Compare sorted versions
        sorted_difference=$(fstrcmp -a $result_sorted $expected_sorted)

        if [[ "$DEBUG" == "FALSE" ]]; then
            # Remove temp files
            rm $result_sorted
            rm $expected_sorted
        else
            echo "# DEBUG: Compare $result_sorted to $expected_sorted"
        fi
        
        # Is result file close enough to expected?
        if [[ "$SHORT_PATH" == "TRUE" ]]; then
            message="Sorted fstrcmp ($sorted_difference) between files $(basename $result_file) and $(basename $expected_file)" 
        else
            message="Sorted fstrcmp ($sorted_difference) between files $result_file and $expected_file" 
        fi
        if (( $( echo "$sorted_difference>$min_score" | bc -l ) )); then
            passed "$message" 2
        else 
            failed "$message" 2
        fi
    fi
}


test_bam_file () {
    module load samtools

    result_file=$1                  # BAM file to compare to expected
    expected_file=$2                # expected BAM file

    samtools view -h $result_file | sort | grep -v '@PG' > $result_file.sorted
    samtools view -h $expected_file | sort | grep -v '@PG' > $expected_file.sorted
    test_text_file $result_file.sorted $expected_file.sorted
}

# Compare Portable Document Format (PDF) file to expected(s)
test_pdf () {
    result_pdf=$1
    expected_pdf=$2
    headless_result=${result_pdf}.headless
    headless_expected=${expected_pdf}.headless
    # First check if the file exists
    if [[ ! -s $result_pdf ]]; then
        failed "File empty or does not exist: $result_pdf"
    else
        tail -n +10 $result_pdf > $headless_result
        tail -n +10 $expected_pdf > $headless_expected

        # See if PDF files after header are identical
        DIFFERENCE=$(get_difference_between $headless_result $headless_expected)

        if [[ "$DEBUG" == "FALSE" ]]; then
            # Clean up temp files
            rm $headless_result $headless_expected
        else
            echo "DEBUG: Compare $headless_result to $headless_expected"
        fi

        # Decide whether we passed, need to check alts, or failed
        if [[ -z "${DIFFERENCE}" ]]; then
            passed_file_is_as_expected "$result_pdf"
        else
            if [[ -s $expected_pdf.alt ]]; then
                test_pdf $result_pdf $expected_pdf.alt
            else
                failed "$result_pdf does not match $expected_pdf"
            fi
        fi
    fi
}

test_value () {
    result=$1
    expected=$2
    message=$3
    if [[ $result == $expected ]]; then
        passed "$message"
    else
        failed "$message (expected: $expected, got: $result)"
    fi
}

approx_csv () {
    original_result=$1
    original_expected=$2
    temp_result=$(mktemp)
    temp_expected=$(mktemp)
    awk 'BEGIN{FS=",";OFS=","} {$2=sprintf("%.2g",$2); $3=sprintf("%.2g",$3); $4=sprintf("%.2g",$4); $5=sprintf("%.2g",$5); $6=sprintf("%.2g",$6); $7=sprintf("%.2g",$7); $8=sprintf("%.2g",$8); $9=sprintf("%.2f",$9); $10=sprintf("%.2f",$10); $11=sprintf("%.2f",$11); $12=sprintf("%.2f",$12); $13=sprintf("%.2g",$13); $14=sprintf("%.2g",$14); $15=sprintf("%.2g",$15); $16=sprintf("%.2g",$16); $17=sprintf("%.2g",$17); $18=sprintf("%.2g",$18); print $0}' $original_result > $temp_result
    awk 'BEGIN{FS=",";OFS=","} {$2=sprintf("%.2g",$2); $3=sprintf("%.2g",$3); $4=sprintf("%.2g",$4); $5=sprintf("%.2g",$5); $6=sprintf("%.2g",$6); $7=sprintf("%.2g",$7); $8=sprintf("%.2g",$8); $9=sprintf("%.2f",$9); $10=sprintf("%.2f",$10); $11=sprintf("%.2f",$11); $12=sprintf("%.2f",$12); $13=sprintf("%.2g",$13); $14=sprintf("%.2g",$14); $15=sprintf("%.2g",$15); $16=sprintf("%.2g",$16); $17=sprintf("%.2g",$17); $18=sprintf("%.2g",$18); print $0}' $original_expected > $temp_expected
    compare_temp_copies $original_result $original_expected $temp_result $temp_expected 

    #Don't fill up TMPDIR
    if [[ "$DEBUG" == "FALSE" ]]; then
        rm -f $temp_result
        rm -f $temp_expected
    else
        echo "# DEBUG: Compare $temp_result to $temp_expected"
    fi
}

test_dir_does_not_exist () {
  dname=$1
  if [ -d $dname ]; then
    failed "$dname exists"
  else
    passed "$dname does not exist"
  fi
}

test_file_exists () {
  filename=$1
  if [ -f $filename ]; then
    passed "$filename exists"
  else
    failed "$filename does not exist"
  fi
}

finished_testing () {
    failed=$(( $TEST_NUM - $TEST_PASSED ))
    if [[ $failed -eq 0 ]]; then
        overall_msg="ALL PASSED"
    else
        overall_msg="FAIL"
    fi
    echo "$overall_msg! (failed/passed/total: $failed/$TEST_PASSED/$TEST_NUM)"

    # Reset test counters
    TEST_NUM=0
    TEST_PASSED=0
}
