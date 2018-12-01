#!/bin/bash -l
#Get the list of variables to analyse
CONFRONTATIONS=$(python2.7 ilamb_setup.py 2>&1)
#
#Option needed to allow ILAMB to run on SPICE
#probably supresses ILAMB trying to access the display
export MPLBACKEND=Agg
#setup ILAMB inputs
#default inputs for optional arguments
echo "###"$LOCAL_INSTALL"###"
echo "###"$REGIONS"###"
echo "###"$CONFRONTATIONS"###"
echo "###"$CLEAN"###"
echo "###"$SKIP_PLOTS"###"
echo "###"$BUILD_DIR"###"
echo "###"$MODELS"###"
echo "###"$MULTI_PROCESSOR"###"

EXTRA=""
if [ $CLEAN == 'true' ]
then
 EXTRA=$EXTRA" --clean"
fi
if [ $SKIP_PLOTS == 'true' ]
then
 EXTRA=$EXTRA" --skip_plots"
fi
if [ ! -z "$REGIONS" ]
then
 EXTRA=$EXTRA" --regions "$REGIONS
fi
if [ ! -z "$CONFRONTATIONS" ]
then
 EXTRA=$EXTRA" --confrontations "$CONFRONTATIONS
fi
if [ ! -z "$BUILD_DIR" ]
then
 EXTRA=$EXTRA" --build_dir "$BUILD_DIR
fi
if [ ! -z "$MODELS" ]
then
 EXTRA=$EXTRA" --models "$MODELS
fi
#save original state so we can get back to it
OLDPATH=$PATH
here=`pwd`
#set up ILAMB python environment
export PATH="/data/users/eroberts/miniconda2/bin:$PATH"
source activate pyEnv_ilamb_2.3_test3
#

if [ $LOCAL_INSTALL == 'true' ]
then
 ILAMB_RUN="$HOME/.local/bin/ilamb-run"
 PYTHONPATH="$PYTHONPATH:/home/h06/eroberts/.local/lib/python2.7/site-packages"
else
 ILAMB_RUN="/data/users/eroberts/miniconda2/envs/pyEnv_ilamb_2.3_test3/bin/ilamb-run"
# PYTHONPATH="$PYTHONPATH:/data/users/eroberts/miniconda2/envs/pyEnv_ilamb_2.3_test3/lib/python2.7/site-packages"
 PATH="/data/users/eroberts/miniconda2/envs/pyEnv_ilamb_2.3_test3/lib/python2.7/site-packages:$PATH"
fi
if [ $MULTI_PROCESSOR == 'true' ]
then
 echo "MULTI_PROCESSOR is TRUE"
 ILAMB_RUN="rose mpi-launch "$ILAMB_RUN
fi
echo "###"$ILAMB_RUN"###"

#
#call ILAMB
export ILAMB_ROOT=$OBS_DIR
echo "$ILAMB_RUN --config $CONFIG_FILE --model_root $MODEL_DIR $EXTRA"
$ILAMB_RUN --config $CONFIG_FILE --model_root $MODEL_DIR $EXTRA
ilambRunExitCode=$?
#important to return to original python environment before exiting
source deactivate
export PATH=$OLDPATH
exit $ilambRunExitCode