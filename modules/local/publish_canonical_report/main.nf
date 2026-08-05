process PUBLISH_CANONICAL_REPORT {
    input:
    path report_files

    output:
    path "report.csv", emit: csv
    path "report.xlsx", emit: xlsx

    script:
    """
    if [ -f nonredundant_report.csv ]; then
        cp nonredundant_report.csv report.csv
        cp nonredundant_report.xlsx report.xlsx
    elif [ -f redundant_report.csv ]; then
        cp redundant_report.csv report.csv
        cp redundant_report.xlsx report.xlsx
    else
        echo "No merged report found" >&2
        exit 1
    fi
    """

    stub:
    """
    touch report.csv
    touch report.xlsx
    """
}
