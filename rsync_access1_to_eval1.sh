rsync -avum \
    rangan@access1.cims.nyu.edu:/data/rangan/dir_bcc/dir_lakcluster_c_dev \
    /home/rangan/dir_bcc/ \
    --include="*/" \
    --include="*.jpg" \
    --include="*.h" \
    --include="*.c" \
    --include="*.in" \
    --include="*.make" \
    --include="*.m" \
    --include="*.sh" \
    --exclude="*" ;
cd /home/rangan/dir_bcc/dir_lakcluster_c_dev/ ;
