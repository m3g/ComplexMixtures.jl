rsync -av --delete ./build/* leandro@leandro:./public_html/m3g/MDDF/
ssh leandro lftp -f /home/leandro/programs/scripts/MDDF_update.lftp
