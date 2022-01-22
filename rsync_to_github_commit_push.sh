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
git add dir_m_dependencies/*.m ;
git add *.make ;
git add *.sh ;
git commit -m "updating dir_lakcluster_c_dev for jeremy elman " ;
git push ;
git pull ;
cd /data/rangan/dir_bcc/dir_lakcluster_c_dev ;

#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
#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
#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
#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 
#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
#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
#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
#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
