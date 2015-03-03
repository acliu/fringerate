#! /bin/bash
for SRC in `ls -d d*` ; do 
    cd $SRC 
    uv2srctrack.py -C psa32_maxred -s $SRC -c 100 *.uvAs
    mv ${SRC}__srctrack.npz ${SRC}_sim.npz
    # why this line?
    rm -rf lst.24559{03.22,46.69}*.uvAs{F,L}
    #uv2srctrack.py -C psa32_maxred -s $SRC -c 100 *.uvAsF
    uv2srctrack.py -C psa32_maxred -s $SRC -c 100 *.uvAsL
    mv ${SRC}__srctrack.npz ${SRC}_fng.npz
    cd .. 
done
