#! /bin/bash
for SRC in `ls -d d*`; do
    echo Working on $SRC
    mdlvis.py -C psa32_maxred -s $SRC lst*.uvA
    fringe_rate_filter_weighted.py -a 0_16 -p I -C psa32_maxred --clean=1 *uvAs
    #plot_uv.py *As{,F} 
    #mkdir d+00
    mv *.uvAs *.uvAsF $SRC
done
