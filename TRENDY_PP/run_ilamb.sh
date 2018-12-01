#!/bin/bash -l
#Get the list of variables to analyse
CONFRONTATIONS=$(python2.7 ilamb_setup.py 2>&1)
#
#Option needed to allow ILAMB to run on SPICE
#probably supresses ILAMB trying to access the display
export MPLBACKEND=Agg
#setup ILAMB inputs
#default inputs for optional arguments
echo "###"$ILAMB_ROOT"###"
echo "###"$REGIONS"###"
echo "###"$CONFRONTATIONS"###"
echo "###"$CLEAN"###"
echo "###"$BUILD_DIR"###"
echo "###"$MODELS"###"
EXTRA=""
if [ $CLEAN == 'true' ]
then
 EXTRA=$EXTRA" --clean"
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
conda_dir=/data/users/eroberts/miniconda2
export PATH="$conda_dir/bin:$PATH"
conda=$conda_dir/bin/conda
source activate env_110917a
#call ILAMB
export ILAMB_ROOT=$ILAMB_ROOT
echo "--config $CONFIG_FILE --model_root $MODEL_DIR $EXTRA"
$ILAMB_ROOT/bin/ilamb-run --config $CONFIG_FILE --model_root $MODEL_DIR $EXTRA
#important to return to original python environment before exiting
source deactivate
export PATH=$OLDPATH