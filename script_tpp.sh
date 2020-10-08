#!/bin/bash

./contact_matrix_asymm DI_ij.dat << EOF
1000
EOF

./tpp contact_corr_MI_2d.xvg contact_map_2gdi_exp_B.xvg > t.xvg 
