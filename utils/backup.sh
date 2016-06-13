#!/bin/bash
tar -cjvf ~/ajm_repo.tar.bz2 /u/fs1/morris/structure_repository
ssh me388@henri.tcm.phy.private.cam.ac.uk 'mv /BIG/me388/ajm_repo.tar.bz2 /BIG/me388/ajm_repo.tar.bz2_bak'
scp ~/ajm_repo.tar.bz2 me388@henri.tcm.phy.private.cam.ac.uk:/BIG/me388
rm ~/ajm_repo.tar.bz2
echo "Backup created successfully."


