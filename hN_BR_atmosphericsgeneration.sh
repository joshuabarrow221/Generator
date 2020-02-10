# Source this script in /genie/app/users/jbarrow/genie-v3/Generator for
# things to build and run correctly on a GENIE gpvm machine
# Source from GENIE 3.0.4 build
source /genie/app/users/jbarrow/genie-v3/genie_env.sh;
rm nohup.out;
cd /pnfs/genie/persistent/users/jbarrow/
# Start atmospheric generation with hN Intranuke and a Bodek-Ritchi relativistic nonlocal Fermi gas model; rotate all events into DUNE coordinates
gevgen_atmo -f HAKKM:/genie/app/users/jbarrow/genie-v3/Generator/hms-ally-20-12-solmax_3FlavOsc.d[12],/genie/app/users/jbarrow/genie-v3/Generator/hms-ally-20-12-solmax_3FlavOsc.d[-12],/genie/app/users/jbarrow/genie-v3/Generator/hms-ally-20-12-solmax_3FlavOsc.d[14],/genie/app/users/jbarrow/genie-v3/Generator/hms-ally-20-12-solmax_3FlavOsc.d[-14],/genie/app/users/jbarrow/genie-v3/Generator/hms-ally-20-12-solmax_3FlavOsc.d[16],/genie/app/users/jbarrow/genie-v3/Generator/hms-ally-20-12-solmax_3FlavOsc.d[-16] -g 1000180400 -n 2000000 -R 0.125237636,-1.57079633,0.0 -E 0.1,100.0 -o NNBarAtm_hN_BR --cross-sections /genie/app/users/jbarrow/genie-v3/Generator/nuall_NNBarAtm_hN_BR.xml --tune J19_00b_00_00b >> NNBarAtm_hN_BR.txt;

# Move all events to proper folder on DUNE gpvms
scp NNBarAtm_hN_BR.100000000.ghep.root jbarrow@dunegpvm08.fnal.gov:/dune/app/users/jbarrow/NEW_WORK/atmospherics/NNBarAtm/hN_BodekRitchie/.;

cd /genie/app/users/jbarrow/genie-v3/Generator/