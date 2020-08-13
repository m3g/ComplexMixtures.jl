rsync -av --delete ./build/* leandro@leandro:./public_html/m3g/ComplexMixtures
ssh leandro lftp -f /home/leandro/programs/scripts/ComplexMixtures_update.lftp
