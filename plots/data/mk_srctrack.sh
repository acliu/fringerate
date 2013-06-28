#! /bin/bash
for SRC in `ls -d d*` ; do 
    cd $SRC 
    #uv2srctrack.py -C psa32_maxred -s $SRC -c 100 *.uvAs
    #mv ${SRC}__srctrack.npz ${SRC}_sim.npz
    rm -rf lst.24559{03.22,46.69}*.uvAsF
    uv2srctrack.py -C psa32_maxred -s $SRC -c 100 *.uvAsF
    mv ${SRC}__srctrack.npz ${SRC}_fng.npz
    cd .. 
done
