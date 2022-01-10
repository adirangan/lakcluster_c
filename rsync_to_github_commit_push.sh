rsync -avum \
    /data/rangan/dir_bcc/dir_lakcluster_c_dev \
    /home/rangan/dir_bcc/ \
    --include="*/" \
    --include="*.h" \
    --include="*.c" \
    --include="*.in" \
    --include="*.make" \
    --include="*.m" \
    --include="*.sh" \
    --exclude="*" ;
cd /data/rangan/dir_bcc/dir_lakcluster_c_dev/ ;
git add dir_h/*.h ;
git add dir_c/*.c ;
git add dir_in/*.in ;
git add dir_m/*.m ;
git add *.make ;
git add *.sh ;
git commit -m "updating dir_lakcluster_c_dev for jeremy elman " ;
git push ;
git pull ;
cd /data/rangan/dir_bcc/dir_lakcluster_c_dev ;
