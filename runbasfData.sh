#!/bin/sh
export USE_GRAND_REPROCESS_DATA=1
export BELLE_MESSAGE_LEVEL=INFO
export LD_LIBRARY_PATH=.:`root-config --libdir`:$LD_LIBRARY_PATH
export BASF_MODULE_DIR=.:$BASF_MODULE_DIR
#setenv BASF_NPROCESS 0

echo    $LD_LIBRARY_PATH
echo    $BASF_MODULE_DIR
echo    $PANTHER_TABLE_DIR


MODULE=handAna
HBKFILE=handana.hbk

basf <<EOF
path create main
path create analysis
path add_module main fix_mdst
path add_condition main >:0:analysis
path add_condition main <=:0:KILL
path add_module analysis ${MODULE}
module put_parameter handAna rfname\dataEx55.root
module put_parameter fix_mdst Make_pi0_option\2
module put_parameter fix_mdst Make_pi0_lower_limit\-5.0
module put_parameter fix_mdst Make_pi0_upper_limit\5.0
module put_parameter fix_mdst Correct_ecl_option\1
initialize
nprocess set 0
histogram define ${HBKFILE}
process_event /pic-disk/projects/belle/ZMDST/DATA/SKIM/HADRON_BJ/EXP55/continuum/HadronBJ-e000055r000793-b20090127_0910.mdst 0
#process_event /pic-disk/projects/belle/ZMDST/DATA/SKIM/HADRON_BJ/EXP55/continuum/HadronBJ-e000055r000794-b20090127_0910.mdst 0
#process_event /pic-disk/projects/belle/ZMDST/DATA/SKIM/HADRON_BJ/EXP55/continuum/HadronBJ-e000055r000795-b20090127_0910.mdst 0
#process_event /pic-disk/projects/belle/ZMDST/DATA/SKIM/HADRON_BJ/EXP55/continuum/HadronBJ-e000055r000796-b20090127_0910.mdst 0
#process_event /pic-disk/projects/belle/ZMDST/DATA/SKIM/HADRON_BJ/EXP55/continuum/HadronBJ-e000055r000797-b20090127_0910.mdst 0
#process_event /pic-disk/projects/belle/ZMDST/DATA/SKIM/HADRON_BJ/EXP55/continuum/HadronBJ-e000055r000798-b20090127_0910.mdst 0

#process_event /pic-disk/projects/belle/ZMDST/DATA/SKIM/HADRON_BJ/EXP55/continuum/HadronBJ-e000055r000799-b20090127_0910.mdst 0
#process_event /pic-disk/projects/belle/ZMDST/DATA/SKIM/HADRON_BJ/EXP55/continuum/HadronBJ-e000055r000800-b20090127_0910.mdst 0
#process_event /pic-disk/projects/belle/ZMDST/DATA/SKIM/HADRON_BJ/EXP55/continuum/HadronBJ-e000055r000801-b20090127_0910.mdst 0
#process_event /pic-disk/projects/belle/ZMDST/DATA/SKIM/HADRON_BJ/EXP55/continuum/HadronBJ-e000055r000802-b20090127_0910.mdst 0
#process_event /pic-disk/projects/belle/ZMDST/DATA/SKIM/HADRON_BJ/EXP55/continuum/HadronBJ-e000055r000805-b20090127_0910.mdst 0
#process_event /pic-disk/projects/belle/ZMDST/DATA/SKIM/HADRON_BJ/EXP55/continuum/HadronBJ-e000055r000806-b20090127_0910.mdst 0
o
utput close
terminate
EOF
