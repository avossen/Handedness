#!/bin/sh
export USE_GRAND_REPROCESS_DATA=1
export BELLE_MESSAGE_LEVEL=INFO
export LD_LIBRARY_PATH=.:`root-config --libdir`:$LD_LIBRARY_PATH
export BASF_Module_Dir=.:$basf_Module_DIR
#setenv BASF_NPROCESS 0

echo    $LD_LIBRARY_PATH
echo    $BASF_MODULE_DIR
echo    $PANTHER_TABLE_DIR


MODULE=handAna
HBKFILE=/dev/null

basf <<EOF
path create main
path create analysis
path add_module main fix_mdst
path add_condition main >:0:analysis
path add_condition main <=:0:KILL
path add_module analysis ${MODULE}
module put_parameter handAna rfname\/people/voss771/mcEx55.root
module put_parameter fix_mdst Make_pi0_option\2
module put_parameter fix_mdst Make_pi0_lower_limit\-20.0
module put_parameter fix_mdst Make_pi0_upper_limit\25.0
module put_parameter fix_mdst Correct_ecl_option\1
initialize
nprocess set 0
histogram define ${HBKFILE}
process_event /pic/projects/belle/ZMDST/DATA/SKIM/HADRON_BJ/EXP55/on_resonance/HadronBJ-e000055r000003-b20090127_0910.mdst 0
#process_event /pic/projects/Belle/ZMDST/MC/MC_4S/MC_4S_EXP55/MC_4S_EXP55/on_resonance/uds/evtgen-uds-00-all-e000055r000508-b20090127_0910.mdst 100000
#process_event /pic/projects/Belle/ZMDST/MC/MC_4S/MC_4S_EXP55/MC_4S_EXP55/on_resonance/uds/evtgen-uds-00-all-e000055r001514-b20090127_0910.mdst-001 10000
#process_event /pic/projects/Belle/ZMDST/MC/MC_4S/MC_4S_EXP55/on_resonance/uds/evtgen-uds-00-all-e000055r000867-b20090127_0910.mdst 50000
#process_event   /pic/projects/Belle/ZMDST/MC/MC_4S/MC_4S_EXP55/on_resonance/uds/evtgen-uds-00-all-e000055r000867-b20090127_0910.mdst 10000
#process_event /pic-disk/projects/belle/ZMDST/MC/MC_4S/MC_4S_EXP43/on_resonance/charm/evtgen-charm-00-all-e000043r000004-b20090127_0910.mdst 10000
#process_event /pic-disk/projects/belle/ZMDST/MC/MC_4S/MC_4S_EXP55/continuum/uds/evtgen-uds-00-all-e000055r001635-b20090127_0910.mdst 100000

#process_event /pic-disk/projects/belle/ZMDST/MC/MC_4S/MC_4S_EXP55/continuum/uds/evtgen-uds-00-all-e000055r001635-b20090127_0910.mdst 0
#process_event /pic-disk/projects/belle/ZMDST/MC/MC_4S/MC_4S_EXP55/continuum/uds/evtgen-uds-00-all-e000055r000801-b20090127_0910.mdst 0
#process_event /pic-disk/projects/belle/ZMDST/MC/MC_4S/MC_4S_EXP55/continuum/uds/evtgen-uds-00-all-e000055r000802-b20090127_0910.mdst
#process_event /pic-disk/projects/belle/ZMDST/MC/MC_4S/MC_4S_EXP55/continuum/uds/evtgen-uds-00-all-e000055r000806-b20090127_0910.mdst
#process_event /pic-disk/projects/belle/ZMDST/MC/MC_4S/MC_4S_EXP55/continuum/uds/evtgen-uds-00-all-e000055r000817-b20090127_0910.mdst
#process_event /pic-disk/projects/belle/ZMDST/MC/MC_4S/MC_4S_EXP55/continuum/uds/evtgen-uds-00-all-e000055r000818-b20090127_0910.mdst
#process_event /pic-disk/projects/belle/ZMDST/MC/MC_4S/MC_4S_EXP55/continuum/uds/evtgen-uds-00-all-e000055r000818-b20090127_0910.mdst-001
#process_event /pic-disk/projects/belle/ZMDST/MC/MC_4S/MC_4S_EXP55/continuum/uds/evtgen-uds-00-all-e000055r000819-b20090127_0910.mdst
#process_event /pic/projects/belle/ZMDST/MC/MC_4S/MC_4S_EXP55/on_resonance/uds/evtgen-uds-00-all-e000055r001429-b20090127_0910.mdst-003 0
#process_event /pic-disk/projects/belle/ZMDST/MC/MC_4S/MC_4S_EXP55/continuum/uds/evtgen-uds-00-all-e000055r000819-b20090127_0910.mdst-001
#process_event /pic-disk/projects/belle/ZMDST/MC/MC_4S/MC_4S_EXP55/continuum/uds/evtgen-uds-00-all-e000055r000820-b20090127_0910.mdst
#/pic/projects/belle/ZMDST/MC/MC_4S/MC_4S_EXP55/on_resonance/uds/evtgen-uds-00-all-e000055r001429-b20090127_0910.mdst-003 
output close
terminate
EOF
