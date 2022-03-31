nice -n 5 matlab -nodisplay -nojvm -nosplash -nodesktop -r \
    "try, run('main.m'), catch, exit(1), end, exit(0);"
echo "matlab exit code: $?"
